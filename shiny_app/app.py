"""
LDT-TargetDB: Ligand-Directed Transfer Radiopharmaceutical Target Database
v2.4 (Major Revision): PU-v6部署分数(44特征: +微环境去卷积6 + P2Rank 2, 54基准AUROC 0.830)
+ 适应症特异性视图 + 失败模式注释 + 口袋3D可视化 + FAQ
"""
import gzip
import json
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px
import streamlit.components.v1 as components

BASE = Path(__file__).resolve().parent.parent
if (BASE / "streamlit_data").exists():
    DATA_DIR = BASE / "streamlit_data"
else:
    DATA_DIR = BASE / "results"

st.set_page_config(page_title="LDT-TargetDB", layout="wide")
st.title("LDT-TargetDB")
st.markdown("**Ligand-Directed Transfer Covalent Radiopharmaceutical Target Prioritization**")
st.caption("Multi-omics + structural chemoproteomics for LDT-NAS radiopharmaceutical target discovery | "
           "[Manuscript preprint](https://github.com/Londger/LDT-TargetDB) | "
           "[FAQ / Help](#faq-help) — see the FAQ tab below")

# ============================================================
# Load Data
# ============================================================
@st.cache_data
def load_data():
    v6_file = DATA_DIR / "enhanced_ranking_v6.csv"
    if v6_file.exists():          # v6: 全量合并表 (PU v6部署分数/微环境/不确定度/HPA蛋白/糖基化/P2Rank/CT.gov)
        return pd.read_csv(v6_file)
    v5_file = DATA_DIR / "enhanced_ranking_v5.csv"
    if v5_file.exists():          # v5回退
        return pd.read_csv(v5_file)
    enhanced = pd.read_csv(DATA_DIR / "enhanced_final_ranking.csv")
    cancer = pd.read_csv(DATA_DIR / "cancer_type_specific_expression.csv")
    covalent = pd.read_csv(DATA_DIR / "covalent_strategy_comparison.csv")
    enhanced = enhanced.merge(cancer[["gene", "best_cancer", "best_log2tpm", "n_cancers_positive"]],
                              on="gene", how="left")
    if "LDT_NAS_score" in covalent.columns and "LDT_NAS_score" not in enhanced.columns:
        enhanced = enhanced.merge(
            covalent[["gene", "LDT_NAS_score", "SuFEx_CTR_score",
                      "Traditional_Acrylamide_score", "best_strategy"]],
            on="gene", how="left")
    v2_file = DATA_DIR / "enhanced_ranking_v2.csv"
    if v2_file.exists():
        v2 = pd.read_csv(v2_file)
        cols = [c for c in ["gene", "enhanced_score_v2", "rank_v2", "eas",
                            "ntrs", "risk_level"] if c in v2.columns]
        enhanced = enhanced.merge(v2[cols], on="gene", how="left")
    return enhanced


@st.cache_data
def load_indication():
    """适应症ITSI宽表 (D10)"""
    p = DATA_DIR / "indication_itsi_wide.csv"
    if not p.exists():
        p = BASE / "results" / "indication_itsi_wide.csv"
    if p.exists():
        return pd.read_csv(p).set_index("gene")
    return None


@st.cache_data
def load_indication_long():
    p = DATA_DIR / "indication_itsi_matrix.csv"
    if not p.exists():
        p = BASE / "results" / "indication_itsi_matrix.csv"
    if p.exists():
        return pd.read_csv(p)
    return None


@st.cache_data
def load_pockets():
    """加载修订版结构分析结果 (top-200胞外口袋)"""
    rev_file = DATA_DIR / "revision" / "revision_structural_analysis.csv"
    if not rev_file.exists():
        rev_file = DATA_DIR / "revision_structural_analysis.csv"
    if rev_file.exists():
        return pd.read_csv(rev_file)
    return None


@st.cache_data
def load_pocket_pdbs():
    """加载预过滤的胞外域PDB文本 (用于3D渲染)"""
    for p in [DATA_DIR / "figures" / "pockets" / "pocket_pdbs.json.gz",
              DATA_DIR / "pocket_pdbs.json.gz"]:
        if p.exists():
            with gzip.open(p, "rt", encoding="utf-8") as f:
                return json.load(f)
    return {}


df = load_data()
pockets_df = load_pockets()
pocket_pdbs = load_pocket_pdbs()
itsi_wide = load_indication()
itsi_long = load_indication_long()

NUC_COLORS = {"LYS": "#2196F3", "CYS": "#FF9800", "TYR": "#4CAF50", "SER": "#9C27B0"}

# ============================================================
# Sidebar
# ============================================================
st.sidebar.header("Filters")
score_col = "enhanced_score_v2" if "enhanced_score_v2" in df.columns else "enhanced_score"
min_score = st.sidebar.slider("Minimum Score", 0.0, 1.0, 0.0, 0.01)
selected_residue = st.sidebar.multiselect(
    "Pocket Nucleophile", ["LYS", "CYS", "TYR", "SER"],
    default=["LYS", "CYS", "TYR", "SER"]
)
strategy_filter = st.sidebar.selectbox(
    "Covalent Strategy", ["All", "LDT_NAS", "SuFEx_CTR", "Traditional_Acrylamide"]
)

# Apply filters
out = df[df[score_col] >= min_score].copy()
if selected_residue:
    out = out[out["top_res_type"].isin(selected_residue)]
if strategy_filter != "All" and "best_strategy" in out.columns:
    out = out[out["best_strategy"] == strategy_filter]

# ============================================================
# Dashboard
# ============================================================
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Targets", len(out))
col2.metric("With Nucleophiles", out["n_nucleophiles"].gt(0).sum())
col3.metric("LYS in Pocket", out["top_res_type"].eq("LYS").sum())
col4.metric("Mean pLDDT", f"{out['mean_plddt'].mean():.0f}" if 'mean_plddt' in out.columns else "N/A")

tab1, tab_ind, tab2, tab3, tab4, tab5 = st.tabs(
    ["Ranking", "Indication view", "Visualization", "Detail", "FAQ / Help", "Download"])

# ============================================================
# Tab: Indication view (D10 - 适应症特异性优先)
# ============================================================
with tab_ind:
    st.subheader("Indication-specific prioritization")
    st.caption(
        "Pan-cancer TSI dilutes indication-specific targets (GPC3-liver rank 1, CEACAM5-pancreatic rank 5, "
        "CLDN18-pancreatic rank 24). "
        "This view ranks targets by the indication-specific expression score (ITSI = log2FC x tumor-positive rate / "
        "sqrt(normal breadth + 1)) computed per cancer type against matched GTEx normal tissues. "
        "Note: CLDN18's gastric ITSI is deeply negative (strong CLDN18 expression in normal gastric mucosa); "
        "indications without TCGA cohorts (e.g. SCLC/DLL3, GEP-NET) cannot be ranked here.")
    if itsi_long is None:
        st.info("Indication matrix not available.")
    else:
        cts = sorted(itsi_long["cancer"].unique())
        sel_ct = st.selectbox("Cancer type (TCGA)", cts,
                              index=cts.index("Pancreatic") if "Pancreatic" in cts else 0)
        sub = itsi_long[itsi_long["cancer"] == sel_ct].sort_values("itsi", ascending=False)
        info_cols = [c for c in ["gene", "pu_rank_v6", "top_res_type", "tier", "label"]
                     if c in df.columns]
        show_i = sub.head(40).merge(df[info_cols], on="gene", how="left")
        if "label" in show_i.columns:
            show_i["known"] = np.where(show_i.label == 1,
                                       "yes" + (f" ({show_i.tier})" if "tier" in show_i else ""), "")
        else:
            show_i["known"] = ""
        out_cols = [c for c in ["gene", "itsi", "fc", "p_tumor", "median_log2", "normal_med",
                                "pu_rank_v6", "top_res_type", "known"] if c in show_i.columns]
        st.dataframe(
            show_i[out_cols],
            use_container_width=True, hide_index=True,
            column_config={
                "gene": "Gene",
                "itsi": st.column_config.NumberColumn("ITSI (indication)", format="%.3f"),
                "fc": st.column_config.NumberColumn("log2FC vs matched normal", format="%.2f"),
                "p_tumor": st.column_config.NumberColumn("Tumor+ rate", format="%.2f"),
                "median_log2": st.column_config.NumberColumn("Median log2TPM", format="%.2f"),
                "normal_med": st.column_config.NumberColumn("Matched normal log2TPM", format="%.2f"),
                "pu_rank_v6": st.column_config.NumberColumn("Deployed rank"),
                "top_res_type": "Best Residue",
                "known": "Known radioligand target",
            })
        st.caption("Note: TCGA does not cover small-cell lung cancer or GEP-NET primaries; "
                   "for such indications, ranked targets reflect the closest available TCGA cohorts and "
                   "supplementary GEO panels (see manuscript Discussion).")

# ============================================================
# Tab 1: Ranking (点击行直接查看详情)
# ============================================================
with tab1:
    st.subheader("Target Ranking")
    st.caption("Click a row to open the gene detail view (or use the Detail tab to search by gene symbol).")

    display_cols = ["gene", "pu_rank_v6", "pu_score_v6", "pu_score_v5",
                    "bag_sd_v6", "pu_score_v5w", "pu_score_v4", "pu_score_v2b",
                    "composite_primary", score_col, "tsi_norm", "ldt_norm",
                    "top_res_type", "n_nucleophiles", "best_cancer", "hpa_tumor_level_max",
                    "tier", "depmap_label", "mean_plddt"]
    display_cols = [c for c in display_cols if c in out.columns]
    show = out[display_cols].head(300)

    event = st.dataframe(
        show, use_container_width=True,
        on_select="rerun", selection_mode="single-row",
        column_config={
            "gene": "Gene",
            "pu_rank_v6": st.column_config.NumberColumn("Rank (deployed)", help="Proteome-wide rank of the deployed v6 score"),
            "pu_score_v6": st.column_config.NumberColumn("PU Score (deployed v6)", format="%.3f", help="Deployed machine-learning prioritization score: PU-HGB ensemble (30 bootstrap bags) on 44 features (36 multi-omics/GEO/HPA/glycosylation + 6 microenvironment deconvolution features [tumor-purity expression attribution, high-purity ITSI, intratumor CV] + 2 P2Rank ligandability features). Out-of-fold on the 54-target benchmark: AUROC 0.830 [0.785-0.868], CV 0.829 ± 0.045, paired vs v5 +0.035 [−0.027,+0.091]; known targets scored out-of-fold"),
            "pu_score_v5": st.column_config.NumberColumn("PU v5 Score", format="%.3f", help="Previous generation: 36 features (27 multi-omics + GEO + HPA protein + glycosylation); 54-target retrain AUROC 0.794 [0.747-0.837]"),
            "bag_sd_v6": st.column_config.NumberColumn("Uncertainty (SD)", format="%.3f", help="Bagged uncertainty: standard deviation of predictions across the 30 bootstrap PU bags of the deployed model. Low SD = stable consensus; high SD = prediction is model-dependent"),
            "pu_score_v5w": st.column_config.NumberColumn("PU v5w Score", format="%.3f", help="Tier-weighted sensitivity variant (v5): positives weighted by clinical evidence tier (approved=3, phase II/III=2, phase I/historical=1); 54-target retrain AUROC 0.792"),
            "pu_score_v4": st.column_config.NumberColumn("PU v4 Score", format="%.3f", help="30-feature generation (27 + GEO indication panels), 47 positives"),
            "pu_score_v2b": st.column_config.NumberColumn("PU v2b Score", format="%.3f", help="10-feature PU ensemble"),
            "composite_primary": st.column_config.NumberColumn("Composite Score", format="%.3f", help="Primary composite score (Methods, Eq. composite)"),
            score_col: st.column_config.NumberColumn("Enhanced Score", format="%.3f", help="Enhanced prioritization score (EAS-integrated)"),
            "tsi_norm": st.column_config.NumberColumn("TSI", format="%.3f"),
            "ldt_norm": st.column_config.NumberColumn("LDT", format="%.3f"),
            "top_res_type": "Best Residue",
            "n_nucleophiles": "# Nuc",
            "best_cancer": "Top Cancer",
            "hpa_tumor_level_max": st.column_config.NumberColumn("HPA protein (max IHC)", format="%.0f", help="Human Protein Atlas tumor IHC: maximum staining level across cancer types (0=not detected, 3=high) - protein-level validation of the RNA signal"),
            "tier": st.column_config.NumberColumn("Evidence tier", help="Clinical evidence tier of known radioligand targets: tier3=approved drug, tier2=phase II/III trials, tier1=phase I/historical"),
            "depmap_label": "DepMap",
            "mean_plddt": st.column_config.NumberColumn("pLDDT", format="%.0f"),
        }, hide_index=True)

    # 点击行 -> 存入session_state, Detail标签页自动显示
    sel = event.selection.get("rows", [])
    if sel:
        sel_gene = str(show.iloc[sel[0]]["gene"])
        st.session_state["selected_gene"] = sel_gene
        st.info(f"Selected: **{sel_gene}** — see the Detail tab for scores, pocket visualization and residue table.")

# ============================================================
# Tab 2: Visualization
# ============================================================
with tab2:
    st.subheader("TSI vs LDT")
    color_map = NUC_COLORS
    fig = px.scatter(out, x="tsi_norm", y="ldt_norm", color="top_res_type",
                     hover_name="gene", title="Tumor Specificity vs LDT Transferability",
                     color_discrete_map=color_map)
    st.plotly_chart(fig, use_container_width=True)

    ca, cb = st.columns(2)
    with ca:
        rc = out["top_res_type"].value_counts()
        fig2 = px.pie(values=rc.values, names=rc.index,
                      title="Nucleophile Distribution", color_discrete_map=color_map)
        st.plotly_chart(fig2, use_container_width=True)
    with cb:
        if "mean_plddt" in out.columns:
            fig3 = px.histogram(out, x="mean_plddt", nbins=30,
                                 title="pLDDT Distribution",
                                 labels={"mean_plddt": "Mean pLDDT"})
            st.plotly_chart(fig3, use_container_width=True)

# ============================================================
# Tab 3: Gene Detail (含3D口袋可视化 + 残基表)
# ============================================================
with tab3:
    st.subheader("Target Detail")
    default_gene = st.session_state.get("selected_gene", "CLDN4")
    search = st.text_input("Gene symbol (or click a row in the Ranking tab)", default_gene).upper().strip()
    hit = out[out["gene"].str.upper() == search]
    if len(hit) > 0:
        r = hit.iloc[0]
        st.markdown(f"### {r['gene']}")
        c1, c2, c3, c4 = st.columns(4)
        deployed = "pu_score_v6" if "pu_score_v6" in r.index else ("pu_score_v5" if "pu_score_v5" in r.index else "pu_score_fused")
        c1.metric("PU Score (deployed)", f"{r.get(deployed, float('nan')):.3f}" if pd.notna(r.get(deployed)) else "N/A")
        c2.metric("Enhanced Score", f"{r.get(score_col, float('nan')):.3f}")
        c3.metric("Composite Score", f"{r['composite_primary']:.3f}" if 'composite_primary' in r.index and pd.notna(r.get('composite_primary')) else "N/A")
        c4.metric("TSI", f"{r['tsi_norm']:.3f}")

        c5, c6, c7, c8 = st.columns(4)
        c5.metric("LDT (norm)", f"{r['ldt_norm']:.3f}")
        if "eas" in r.index and pd.notna(r.get("eas")):
            c6.metric("EAS", f"{r['eas']:.1f}")
        if "risk_level" in r.index and pd.notna(r.get("risk_level")):
            c7.metric("Normal Tissue Risk", str(r["risk_level"]))
        if "mean_plddt" in r.index and pd.notna(r.get("mean_plddt")):
            c8.metric("Mean pLDDT", f"{r['mean_plddt']:.0f}")

        st.markdown(f"**Nucleophile:** {r.get('top_res_type', '?')} | "
                    f"**Top Cancer:** {r.get('best_cancer', '?')} ({r.get('best_log2tpm', '?')} log2TPM) | "
                    f"**DepMap:** {r.get('depmap_label', '?')} (Chronos {r.get('depmap_score', float('nan')):.3g})"
                    + (f" | **NTRS:** {r['ntrs']:.2f}" if "ntrs" in r.index and pd.notna(r.get("ntrs")) else "")
                    + (" | :star: Known radioligand target" if r.get("known_target") == 1 else ""))

        # ---- v5/v6: 不确定度 + 蛋白层/糖基化/微环境证据 + 临床转化 + 失败模式 (D11) ----
        ev1, ev2, ev3, ev4, ev5, ev6, ev7 = st.columns(7)
        if "bag_sd_v6" in r.index and pd.notna(r.get("bag_sd_v6")):
            ev1.metric("Prediction uncertainty (bag SD)",
                       f"{r['bag_sd_v6']:.3f}",
                       help="Standard deviation of the score across the 30 bootstrap PU bags of the deployed v6 model. "
                            "Known targets typically show SD ~0.04; values >0.15 indicate the rank is model-dependent.")
        elif "bag_sd_v5" in r.index and pd.notna(r.get("bag_sd_v5")):
            ev1.metric("Prediction uncertainty (bag SD)",
                       f"{r['bag_sd_v5']:.3f}",
                       help="Standard deviation of the score across the 30 bootstrap PU bags. "
                            "Low SD = stable consensus; values >0.15 indicate the rank is model-dependent.")
        if "hpa_tumor_level_max" in r.index and pd.notna(r.get("hpa_tumor_level_max")):
            lvl_names = {0.0: "Not detected", 1.0: "Low", 2.0: "Medium", 3.0: "High"}
            ev2.metric("HPA tumor protein (max IHC)",
                       lvl_names.get(float(r["hpa_tumor_level_max"]), r["hpa_tumor_level_max"]),
                       help=f"Protein-level validation: maximum tumor IHC staining across {int(r.get('n_cancers_med_high', 0))} cancer types with medium/high staining. Source: Human Protein Atlas v22 pathology.")
        if "n_glyco_EC" in r.index and pd.notna(r.get("n_glyco_EC")):
            ev3.metric("EC N-glycosylation sites",
                       int(r["n_glyco_EC"]),
                       help="UniProt-annotated extracellular N-linked glycosylation sites. Densely glycosylated pockets may shield residues from both antibody fragments and small-molecule ligands (glycan-near-pocket ratio in the download file).")
        if "ctgov_n_trials" in r.index and pd.notna(r.get("ctgov_n_trials")):
            ev4.metric("CT.gov radioligand trials",
                       int(r["ctgov_n_trials"]),
                       help=f"Highest interventional phase: {r.get('ctgov_max_phase', 'N/A')}. Source: ClinicalTrials.gov API v2 (see Supplementary Table S6).")
        if "p2rank_top_prob" in r.index and pd.notna(r.get("p2rank_top_prob")):
            ev5.metric("P2Rank ligandability",
                       f"{r['p2rank_top_prob']:.3f}",
                       help="Machine-learning pocket prediction (P2Rank 2.5.1) on the full AlphaFold structure: probability that the top-scoring pocket is ligandable. Provided as structural annotation for binder/warhead design; note that classical small-molecule ligandability does not discriminate known radioligand targets, which often bind via extended surfaces.")
        if "p2rank_filt_top_prob" in r.index and pd.notna(r.get("p2rank_filt_top_prob")):
            ev6.metric("P2Rank (EC-filtered)",
                       f"{r['p2rank_filt_top_prob']:.3f}",
                       help=f"P2Rank on the refined structure (pLDDT<50, signal peptide and non-extracellular residues removed; {int(r.get('p2rank_filt_n_pockets', 0))} pockets with probability>=0.5). Consistent with the refined pocket analysis used for ranking.")
        if "purity_corr_max" in r.index and pd.notna(r.get("purity_corr_max")):
            ev7.metric("Tumor-purity corr (max)",
                       f"{r['purity_corr_max']:.2f}",
                       help="Microenvironment deconvolution feature of the deployed v6 model: maximum Spearman correlation between expression and a tumor-purity proxy across TCGA cancer types. Positive = expression attributable to tumor cells rather than infiltrating immune/stromal cells (favorable for radioligand specificity); negative or near-zero = expression may reflect microenvironment contamination.")

        if "tier" in r.index and pd.notna(r.get("tier")):
            tier_txt = {"tier3": "tier 3 - approved radiopharmaceutical",
                        "tier2": "tier 2 - phase II/III trials",
                        "tier1": "tier 1 - phase I / historical"}.get(str(r["tier"]), str(r["tier"]))
            st.markdown(f":star: **Known radioligand target ({tier_txt})**")
        if "failure_note" in r.index and pd.notna(r.get("failure_note")):
            st.warning(f"**Known failure mode:** this target is under-ranked by bulk-expression features - "
                       f"{r['failure_note']}. Consider the indication view and GEO/single-cell evidence.")
        if "subsystem" in r.index and pd.notna(r.get("subsystem")) and r.get("subsystem") == "blood":
            st.caption("Note: with the deployed PU-v6 model, blood-disseminated targets are ranked "
                       "comparably to solid-tumor targets (subset AUROC 0.846 vs 0.827); circulating-cell "
                       "indications without TCGA cohorts remain annotated as failure modes.")

        # ---- 3D口袋可视化 (预过滤胞外域PDB + 口袋残基高亮) ----
        if pockets_df is not None and len(pocket_pdbs) > 0:
            prow = pockets_df[pockets_df["gene"].str.upper() == search]
            if len(prow) > 0:
                r_pdb = pocket_pdbs.get(prow.iloc[0]["gene"])
            else:
                r_pdb = None
            if len(prow) > 0 and r_pdb:
                pr = prow.iloc[0]
                try:
                    pocket_res = json.loads(pr["pocket_residues"]) if pd.notna(pr.get("pocket_residues")) else []
                except Exception:
                    pocket_res = []
                if pocket_res and r_pdb:
                    st.markdown("#### Extracellular pocket (3D)")
                    st.caption("Cartoon: extracellular high-confidence region (pLDDT ≥ 50). "
                               "Gray sticks: pocket-lining residues. Colored sticks: LDT-qualifying "
                               "nucleophiles (Lys blue, Cys orange, Tyr green, Ser purple). "
                               "Drag to rotate; scroll to zoom.")
                    try:
                        import py3Dmol
                        pocket_nums = set(p["res_num"] for p in pocket_res)
                        view = py3Dmol.view(width=700, height=480)
                        view.addModel(r_pdb, "pdb")
                        view.setStyle({"cartoon": {"color": "lightblue", "opacity": 0.5}})
                        view.addStyle({"resi": list(pocket_nums)},
                                      {"stick": {"color": "gray", "radius": 0.15}})
                        for p in pocket_res:
                            if p["res_type"] in NUC_COLORS:
                                view.addStyle({"resi": p["res_num"]},
                                              {"stick": {"color": NUC_COLORS[p["res_type"]], "radius": 0.3}})
                        view.zoomTo()
                        components.html(view._make_html(), height=500, scrolling=False)
                    except ImportError:
                        st.warning("py3Dmol not available; showing pocket residues in the table below.")

                    # ---- 残基表 (含pKa/SASA, 来自top-3口袋明细) ----
                    st.markdown("#### Pocket residues")
                    try:
                        top3 = json.loads(pr["top3_pockets"]) if pd.notna(pr.get("top3_pockets")) else []
                    except Exception:
                        top3 = []
                    rows = []
                    for pk in top3:
                        for nuc in pk.get("nucleophiles", []):
                            rows.append({
                                "Pocket rank": pk["rank"],
                                "Pocket volume (Å³)": pk.get("volume"),
                                "Residue": nuc.get("res", ""),
                                "Type": nuc.get("type", ""),
                                "pKa": nuc.get("pka"),
                                "Rel. SASA (%)": nuc.get("sasa_rel"),
                                "LDT residue score": nuc.get("score"),
                            })
                    if rows:
                        st.dataframe(pd.DataFrame(rows).sort_values(
                            ["Pocket rank", "LDT residue score"], ascending=[True, False]),
                            use_container_width=True, hide_index=True)
                    else:
                        st.caption(f"Largest pocket: {int(len(pocket_res))} residues "
                                   f"(volume {pr.get('max_pocket_volume', '?')} Å³, "
                                   f"mean pocket pLDDT {pr.get('mean_pocket_plddt', '?')}). "
                                   "No LDT-qualifying nucleophile in the top-3 pockets.")
                    st.caption(f"Pocket pLDDT: {pr.get('mean_pocket_plddt', 'N/A')} | "
                               f"Best nucleophile: {pr.get('top_nucleophile', 'N/A')} | "
                               f"pKa {pr.get('top_pka', 'N/A')} | Rel. SASA {pr.get('top_sasa', 'N/A')}%")
            elif len(prow) > 0:
                st.info("No detectable extracellular pocket for this target after filtering "
                        "(pLDDT ≥ 50, extracellular residues only).")

        # AlphaFold直链 (按accession, 不按名称搜索)
        uniprot = r.get('uniprot_id', '')
        if pd.notna(uniprot) and uniprot:
            st.markdown(f"[AlphaFold Structure ({uniprot})](https://alphafold.ebi.ac.uk/entry/{uniprot}) | "
                        f"[UniProt](https://www.uniprot.org/uniprotkb/{uniprot})")

        # Cancer-type specific expression
        ct_data = pd.read_csv(DATA_DIR / "cancer_type_specific_expression.csv")
        ct_hit = ct_data[ct_data["gene"].str.upper() == search]
        if len(ct_hit) > 0:
            st.markdown(f"**Best Cancer Types:** {ct_hit['best_cancer'].values[0]} ({ct_hit['best_log2tpm'].values[0]}), "
                        f"{ct_hit['second_cancer'].values[0]} ({ct_hit['second_log2tpm'].values[0]}), "
                        f"{ct_hit['third_cancer'].values[0]} ({ct_hit['third_log2tpm'].values[0]})")
    else:
        st.warning(f"'{search}' not found.")

# ============================================================
# Tab 4: FAQ / Help (审稿人要求: 缩写与评分说明)
# ============================================================
with tab4:
    st.subheader("FAQ / Help")
    st.markdown(
        """
#### Abbreviations

| Term | Meaning |
|---|---|
| **LDT** | Ligand-Directed Transfer — traceless covalent labeling of endogenous proteins via N-acyl-N-alkyl/aryl sulfonamide (NASA/ArNASA) warheads |
| **TSI** | Tumor Specificity Index — tumor-selective expression score: log2 fold change × tumor positive rate / √(GTEx normal samples above threshold + 1) |
| **LTS** | LDT Transferability Score — best nucleophile score within detected pockets: max(w × pₛ × sₛ × b) |
| **EAS** | Extracellular Accessibility Score — confidence of extracellular localization (1.0 = UniProt-confirmed extracellular topological domain; 0.5 = topology without confirmed EC domain, or neutral default where topology was not systematically assessed) |
| **NTRS** | Normal Tissue Risk Score — maximum organ-weighted HPA staining level across seven radiosensitive organs (higher = riskier) |
| **V_pocket** | Maximum detected pocket volume, min–max normalized to [0, 1] |
| **pLDDT** | Predicted local distance difference test — per-residue AlphaFold confidence (0–100) |
| **TPM** | Transcripts Per Million — RNA-seq normalization unit |
| **DepMap Chronos** | CRISPR gene essentiality score (more negative = more essential) |

#### Scores

- **Primary composite score** (proteome-wide ranking):
  `0.35 × TSI + 0.30 × V_pocket + 0.25 × LTS + 0.10 × (1 − pLDDT_below70/100)`
- **Enhanced prioritization score** (final target selection):
  `0.25 × TSI + 0.20 × V_pocket + 0.20 × LTS × EAS_penalty + 0.10 × pLDDT_term + 0.05 × DepMap + 0.20 × EAS`
- **PU score (deployed, v6)** — the primary machine-learning prioritization score: a PU-HGB
  ensemble (30 bootstrap bags) on 44 features. Building on the v5 feature set (27 multi-omics
  + 3 GEO indication panels + 3 HPA tumor-IHC protein levels + 3 UniProt N-glycosylation
  site counts), v6 adds 6 microenvironment deconvolution features (tumor-purity expression
  attribution, high-purity-subsample ITSI, intratumor expression CV) and 2 P2Rank ML
  ligandability features. Trained on 54 known radioligand targets versus unlabeled surface
  proteins; known targets are scored strictly out-of-fold. Out-of-fold AUROC
  0.830 [0.785–0.868] (CV 0.829 ± 0.045; paired ΔAUROC vs v5 +0.035
  [−0.027, +0.091], P(v6 > v5) = 0.864; on the 47-target pre-expansion benchmark
  +0.032 [+0.009, +0.056], P = 0.997). The microenvironment features correct the bulk-RNA
  dilution of tumor-cell-specific targets (e.g. SSTR2 rank 1,347 → 646, GPC3 → 250,
  CLDN18 → 222 on the 54-target benchmark).
- **Prediction uncertainty (bag SD)** — standard deviation of the score across the 30 bootstrap
  PU bags of the deployed model. Known targets show SD ~0.04 vs ~0.19 for unlabeled proteins;
  a high SD flags candidates whose rank is model-dependent (use the full score distribution,
  not the point rank).
- **Evidence tier** — clinical-trial annotation of the 54 known targets
  (tier 3: n = 5 approved; tier 2: n = 34; tier 1: n = 15):
  tier 3 = approved radiopharmaceutical (e.g. PSMA, SSTR2, CD20), tier 2 = phase II/III trials,
  tier 1 = phase I / historical; see Supplementary Table S6 for the full trial list.
- **Indication view** — per-cancer-type ITSI ranking (32 TCGA types vs matched GTEx normals).
  It recovers indication-specific targets that the pan-cancer TSI dilutes (e.g. GPC3-liver,
  CLDN18-pancreatic, CEACAM5-pancreatic). TCGA does not cover SCLC or GEP-NET primaries;
  these blind spots are annotated on the detail page.
- **Known failure modes** — (1) circulating-cell blood indications without TCGA cohorts remain
  harder to represent (blood-lineage targets overall rank comparably to solid-tumor targets
  under PU-v6: subset AUROC 0.846 vs 0.827); (2) targets whose indications lack TCGA cohorts
  (GEP-NET, insulinoma, SCLC) rely on GEO supplementation, microenvironment deconvolution, and
  single-cell evidence; (3) glycan-shielded pockets may hide ligandable nucleophiles (see EC
  N-glycosylation metric).
- **PU v5 score**: previous generation (36 features; scores shown are the 47-target-era
  retraining, AUROC 0.783 [0.727-0.835]; the final 54-target retrain reaches
  AUROC 0.794 [0.747-0.837]); its
  tier-weighted variant (v5w, 54-target retrain AUROC 0.792) is retained for reference.
- **P2Rank ligandability**: probability that the top-scoring pocket is a small-molecule-binding
  site (P2Rank 2.5.1, random-forest consensus; full AlphaFold structure, plus an
  extracellular-filtered variant for top-200 targets). Provided as structural annotation for
  binder/warhead design. Note: P2Rank ligandability does NOT discriminate known radioligand
  targets (AUROC 0.43 as a standalone ranker) — radioligand targets typically bind ligands
  via extended extracellular surfaces rather than classical druggable pockets, which is why
  the deployed score uses nucleophile-centric LDT chemistry instead.
- **PU v2 score** (positive-unlabeled ensemble): bagged logistic regression + gradient boosting
  integrating indication-aware expression features (per-cancer-type tumor-vs-matched-normal maxITSI),
  pan-cancer TSI, pocket volume, LTS and pLDDT. Known radioligand targets are scored out-of-fold
  (no target is scored by a model trained on itself).
- **LDT residue score** = `w × pₛ × sₛ × b`, where w is the residue type weight
  (Lys 1.0, Cys 0.3, Tyr 0.2, Ser 0.1), pₛ = exp(−|pKa − 7.4|/2) is the pKa reactivity score,
  sₛ = min(relative SASA / 50%, 1) is the solvent accessibility score, and b = 2 for Lys
  (proximity boost).

#### Pocket detection

Pockets were detected with pyKVFinder on AlphaFold structures filtered to extracellular,
high-confidence residues only (pLDDT ≥ 50, signal peptides and non-extracellular regions removed).
For the top 200 targets, the three largest pockets per protein are scored.

#### Two-stage workflow

1. **Stage 1 — expression filter:** TSI ranking identifies tumor-associated surface proteins.
2. **Stage 2 — structural/chemical refinement:** pocket geometry and LDT chemistry
   compatibility help prioritize among expression-equivalent candidates.

#### Data & citation

All scores, the pocket residue composition and pre-filtered extracellular structures for the
top 200 targets can be downloaded from the Download tab. If you use LDT-TargetDB in your
research, please cite the LDT-TargetDB manuscript (Bioinformatics Advances, BIOADV-2026-446).
        """
    )

# ============================================================
# Tab 5: Download
# ============================================================
with tab5:
    st.subheader("Download")
    csv = out.to_csv(index=False)
    st.download_button("Download filtered ranking (CSV)", csv, "LDT_TargetDB.csv", "text/csv")
    if pockets_df is not None:
        slim = pockets_df[[c for c in ["gene", "n_pockets", "max_pocket_volume", "n_nucleophiles",
                                        "ldt_score", "ldt_score_largest", "top_res_type",
                                        "best_pocket_rank", "best_pocket_volume", "ligandable_ldt",
                                        "mean_pocket_plddt", "mean_plddt_kept", "top_nucleophile",
                                        "top_pka", "top_sasa"]
                           if c in pockets_df.columns]]
        st.download_button("Download top-200 extracellular pocket analysis (CSV)",
                           slim.to_csv(index=False), "LDT_TargetDB_pockets.csv", "text/csv")
    if itsi_wide is not None:
        st.download_button("Download per-cancer-type ITSI matrix (CSV)",
                           itsi_wide.to_csv(), "LDT_TargetDB_ITSI_matrix.csv", "text/csv")

st.sidebar.markdown("---")
st.sidebar.markdown("**LDT-TargetDB v2.4** (Major Revision)")
st.sidebar.markdown("New: PU-v6 deployed score (44 features incl. microenvironment deconvolution "
                    "and P2Rank ligandability, 54-benchmark AUROC 0.830), bagged uncertainty, "
                    "indication-specific view, failure-mode annotations, CT.gov trial counts")
st.sidebar.markdown("[GitHub](https://github.com/Londger/LDT-TargetDB)")
