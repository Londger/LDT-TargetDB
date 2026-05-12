"""
LDT-TargetDB: Ligand-Directed Transfer Radiopharmaceutical Target Database
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px

BASE = Path(__file__).resolve().parent.parent
if (BASE / "streamlit_data").exists():
    DATA_DIR = BASE / "streamlit_data"
else:
    DATA_DIR = BASE / "results"

st.set_page_config(page_title="LDT-TargetDB", layout="wide")
st.title("LDT-TargetDB")
st.markdown("**Ligand-Directed Transfer Covalent Radiopharmaceutical Target Prioritization**")
st.caption("Multi-omics + structural chemoproteomics for LDT-NAS radiopharmaceutical target discovery")

# ============================================================
# Load Data
# ============================================================
@st.cache_data
def load_data():
    enhanced = pd.read_csv(DATA_DIR / "enhanced_final_ranking.csv")
    cancer = pd.read_csv(DATA_DIR / "cancer_type_specific_expression.csv")
    covalent = pd.read_csv(DATA_DIR / "covalent_strategy_comparison.csv")
    enhanced = enhanced.merge(cancer[["gene","best_cancer","best_log2tpm","n_cancers_positive"]],
                              on="gene", how="left")
    if "LDT_NAS_score" in covalent.columns and "LDT_NAS_score" not in enhanced.columns:
        enhanced = enhanced.merge(
            covalent[["gene","LDT_NAS_score","SuFEx_CTR_score",
                      "Traditional_Acrylamide_score","best_strategy"]],
            on="gene", how="left")
    return enhanced

df = load_data()

# ============================================================
# Sidebar
# ============================================================
st.sidebar.header("Filters")
min_score = st.sidebar.slider("Minimum Score", 0.0, 1.0, 0.3, 0.05)
selected_residue = st.sidebar.multiselect(
    "Pocket Nucleophile", ["LYS", "CYS", "TYR", "SER"],
    default=["LYS", "CYS", "TYR", "SER"]
)
strategy_filter = st.sidebar.selectbox(
    "Covalent Strategy", ["All", "LDT_NAS", "SuFEx_CTR", "Traditional_Acrylamide"]
)

# Apply filters
out = df[df["enhanced_score"] >= min_score]
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

tab1, tab2, tab3, tab4 = st.tabs(["Ranking", "Visualization", "Detail", "Download"])

with tab1:
    st.subheader("Target Ranking")
    display_cols = ["gene", "enhanced_score", "tsi_norm", "ldt_norm",
                     "top_res_type", "n_nucleophiles", "best_cancer",
                     "depmap_label", "mean_plddt"]
    display_cols = [c for c in display_cols if c in out.columns]
    st.dataframe(
        out[display_cols].head(100), use_container_width=True,
        column_config={
            "gene": "Gene",
            "enhanced_score": st.column_config.NumberColumn("Score", format="%.3f"),
            "tsi_norm": st.column_config.NumberColumn("TSI", format="%.3f"),
            "ldt_norm": st.column_config.NumberColumn("LDT", format="%.3f"),
            "top_res_type": "Best Residue",
            "n_nucleophiles": "# Nuc",
            "best_cancer": "Top Cancer",
            "depmap_label": "DepMap",
            "mean_plddt": st.column_config.NumberColumn("pLDDT", format="%.0f"),
        }, hide_index=True)

with tab2:
    st.subheader("TSI vs LDT")
    color_map = {"LYS": "#e74c3c", "CYS": "#f39c12", "TYR": "#3498db", "SER": "#95a5a6"}
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

with tab3:
    st.subheader("Target Detail")
    search = st.text_input("Gene symbol", "CLDN4").upper()
    hit = out[out["gene"].str.upper() == search]
    if len(hit) > 0:
        r = hit.iloc[0]
        st.markdown(f"### {r['gene']}")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Score", f"{r['enhanced_score']:.3f}")
        c2.metric("TSI", f"{r['tsi_norm']:.3f}")
        c3.metric("LDT", f"{r['ldt_norm']:.3f}")
        c4.metric("Pockets", int(r.get('n_pockets', 0)))
        st.markdown(f"**Nucleophile:** {r.get('top_res_type','?')}")
        st.markdown(f"**Top Cancer:** {r.get('best_cancer','?')} ({r.get('best_log2tpm','?')} log2TPM)")
        st.markdown(f"**DepMap:** {r.get('depmap_label','?')} ({r.get('depmap_score','?')})")
        st.markdown(f"**pLDDT:** {r.get('mean_plddt','?')}")
        uniprot = r.get('uniprot_id', '')
        if pd.notna(uniprot) and uniprot:
            st.markdown(f"[AlphaFold Structure](https://alphafold.ebi.ac.uk/entry/{uniprot})")
        # Cancer-type specific expression
        ct_data = pd.read_csv(DATA_DIR / "cancer_type_specific_expression.csv")
        ct_hit = ct_data[ct_data["gene"].str.upper() == search]
        if len(ct_hit) > 0:
            st.markdown(f"**Best Cancer Types:** {ct_hit['best_cancer'].values[0]} ({ct_hit['best_log2tpm'].values[0]}), "
                        f"{ct_hit['second_cancer'].values[0]} ({ct_hit['second_log2tpm'].values[0]}), "
                        f"{ct_hit['third_cancer'].values[0]} ({ct_hit['third_log2tpm'].values[0]})")
    else:
        st.warning(f"'{search}' not found.")

with tab4:
    st.subheader("Download")
    csv = out.to_csv(index=False)
    st.download_button("Download CSV", csv, "LDT_TargetDB.csv", "text/csv")

st.sidebar.markdown("---")
st.sidebar.markdown("**LDT-TargetDB v1.0**")
st.sidebar.markdown("[GitHub](https://github.com/Londger/LDT-TargetDB)")
