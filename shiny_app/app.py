"""
LDT-TargetDB: Ligand-Directed Transfer Radiopharmaceutical Target Database
Streamlit Web Application
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

BASE = Path(__file__).resolve().parent.parent
# Use streamlit_data/ if available (deployment), else results/
if (BASE / "streamlit_data").exists():
    DATA_DIR = BASE / "streamlit_data"
else:
    DATA_DIR = BASE / "results"

st.set_page_config(page_title="LDT-TargetDB", layout="wide")
st.title("🔬 LDT-TargetDB")
st.markdown("**Ligand-Directed Transfer Covalent Radiopharmaceutical Target Prioritization**")
st.caption("Integrated multi-omics + structural chemoproteomics for LDT-NAS radiopharmaceutical target discovery")

# ============================================================
# Load Data
# ============================================================
@st.cache_data
def load_data():
    enhanced = pd.read_csv(DATA_DIR / "enhanced_final_ranking.csv")
    covalent = pd.read_csv(DATA_DIR / "covalent_strategy_comparison.csv")
    # Merge covalent if not already in enhanced
    if "LDT_NAS_score" in covalent.columns and "LDT_NAS_score" not in enhanced.columns:
        enhanced = enhanced.merge(
            covalent[["gene", "LDT_NAS_score", "SuFEx_CTR_score",
                       "Traditional_Acrylamide_score", "best_strategy"]],
            on="gene", how="left"
        )
    return enhanced

df = load_data()

# ============================================================
# Sidebar Filters
# ============================================================
st.sidebar.header("Filters")
min_score = st.sidebar.slider("Minimum Enhanced Score", 0.0, 1.0, 0.3, 0.05)
min_tsi = st.sidebar.slider("Minimum TSI", 0.0, 1.0, 0.0, 0.05)
selected_residue = st.sidebar.multiselect(
    "LDT Target Residue", ["LYS", "TYR", "HIS", "SER"], default=["LYS", "TYR", "HIS", "SER"]
)
strategy_filter = st.sidebar.selectbox(
    "Covalent Strategy", ["All", "LDT_NAS", "SuFEx_CTR", "Traditional_Acrylamide"]
)

# Apply filters
filtered = df[df["enhanced_score"] >= min_score]
if min_tsi > 0:
    filtered = filtered[filtered["tsi_norm"] >= min_tsi]
if selected_residue:
    filtered = filtered[filtered["top_res_type"].isin(selected_residue)]
if strategy_filter != "All" and "best_strategy" in filtered.columns:
    filtered = filtered[filtered["best_strategy"] == strategy_filter]

# ============================================================
# Main Dashboard
# ============================================================
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Targets", len(filtered))
col2.metric("With LDT Nucleophiles", filtered["n_nucleophiles"].gt(0).sum())
col3.metric("LYS in Pocket", filtered["top_res_type"].eq("LYS").sum())
col4.metric("Mean pLDDT", f"{filtered['mean_plddt'].mean():.0f}" if 'mean_plddt' in filtered.columns else "N/A")

# Tabs
tab1, tab2, tab3, tab4 = st.tabs(["📊 Ranking", "📈 Visualization", "🔍 Target Detail", "📥 Download"])

with tab1:
    st.subheader("Target Ranking")
    display_cols = ["gene", "enhanced_score", "tsi_norm", "ldt_norm",
                     "top_res_type", "n_nucleophiles", "top1_cancer", "depmap_label", "mean_plddt"]
    display_cols = [c for c in display_cols if c in filtered.columns]
    st.dataframe(
        filtered[display_cols].head(100),
        use_container_width=True,
        column_config={
            "gene": "Gene",
            "enhanced_score": st.column_config.NumberColumn("Score", format="%.3f"),
            "tsi_norm": st.column_config.NumberColumn("TSI", format="%.3f"),
            "ldt_norm": st.column_config.NumberColumn("LDT", format="%.3f"),
            "top_res_type": "Best Residue",
            "n_nucleophiles": "# Nucleophiles",
            "top1_cancer": "Top Cancer Type",
            "depmap_label": "DepMap",
            "mean_plddt": st.column_config.NumberColumn("pLDDT", format="%.0f"),
        },
        hide_index=True,
    )

with tab2:
    st.subheader("Score Distribution")

    # Scatter: TSI vs LDT
    fig = px.scatter(
        filtered, x="tsi_norm", y="ldt_norm",
        color="top_res_type", hover_name="gene",
        title="TSI vs LDT Transferability",
        labels={"tsi_norm": "Tumor Specificity (TSI)", "ldt_norm": "LDT Score"},
        color_discrete_map={"LYS": "#e74c3c", "TYR": "#3498db", "HIS": "#2ecc71", "SER": "#95a5a6"}
    )
    st.plotly_chart(fig, use_container_width=True)

    # Residue type distribution
    col_a, col_b = st.columns(2)
    with col_a:
        res_counts = filtered["top_res_type"].value_counts()
        fig2 = px.pie(values=res_counts.values, names=res_counts.index,
                       title="Nucleophile Distribution in Best Pocket")
        st.plotly_chart(fig2, use_container_width=True)

    with col_b:
        if "best_strategy" in filtered.columns:
            strat_counts = filtered["best_strategy"].value_counts()
            fig3 = px.bar(x=strat_counts.index, y=strat_counts.values,
                           title="Optimal Covalent Strategy per Target",
                           labels={"x": "Strategy", "y": "# Targets"})
            st.plotly_chart(fig3, use_container_width=True)

    # pLDDT distribution
    if "mean_plddt" in filtered.columns:
        fig4 = px.histogram(filtered, x="mean_plddt", nbins=30,
                             title="pLDDT Confidence Distribution",
                             labels={"mean_plddt": "Mean pLDDT"})
        st.plotly_chart(fig4, use_container_width=True)

with tab3:
    st.subheader("Target Detail Lookup")
    search_gene = st.text_input("Enter gene symbol", "CLDN4").upper()
    target = filtered[filtered["gene"].str.upper() == search_gene]
    if len(target) > 0:
        row = target.iloc[0]
        st.markdown(f"### {row['gene']}")
        st.markdown(f"**{row.get('description','')}**")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Enhanced Score", f"{row['enhanced_score']:.3f}")
        c2.metric("TSI", f"{row['tsi_norm']:.3f}")
        c3.metric("LDT Score", f"{row['ldt_norm']:.3f}")
        c4.metric("Pockets", int(row.get('n_pockets', 0)))

        st.markdown(f"**Best Nucleophile:** {row['top_res_type']} (pKa={row.get('top_res_pka','N/A')})")
        st.markdown(f"**Top Cancer Type:** {row.get('top1_cancer','N/A')}")
        st.markdown(f"**DepMap Status:** {row.get('depmap_label','N/A')} (score={row.get('depmap_score','N/A')})")
        st.markdown(f"**pLDDT:** {row.get('mean_plddt','N/A')} (±{row.get('plddt_below_70','N/A')}% <70)")

        # AlphaFold structure link
        uniprot = row.get('uniprot_id', '')
        if pd.notna(uniprot) and uniprot:
            st.markdown(f"[View AlphaFold Structure](https://alphafold.ebi.ac.uk/entry/{uniprot})")
    else:
        st.warning(f"Gene '{search_gene}' not found in filtered results.")

with tab4:
    st.subheader("Download Data")
    csv = filtered.to_csv(index=False)
    st.download_button("Download CSV", csv, "LDT_TargetDB_results.csv", "text/csv")
    st.caption("Full dataset including all scoring dimensions.")

st.sidebar.markdown("---")
st.sidebar.markdown("**LDT-TargetDB v1.0**")
st.sidebar.markdown("Developed for covalent radiopharmaceutical target discovery")
st.sidebar.markdown("[GitHub Repository](#)")
