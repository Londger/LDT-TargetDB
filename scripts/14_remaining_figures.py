"""
补全Fig6-8 + Table1-4
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
import json

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
FIG_DIR = BASE / "results" / "figures"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams.update({'font.family': 'Arial', 'font.size': 10,
                     'axes.titlesize': 12, 'axes.labelsize': 11})

# Load data
enhanced = pd.read_csv(BASE / "results/enhanced_final_ranking.csv")
ct_df = pd.read_csv(BASE / "results/cancer_type_specific_expression.csv")
tool_comp = pd.read_csv(BASE / "results/tool_comparison.csv")
known = pd.read_csv(BASE / "results/known_targets_validation.csv")
covalent = pd.read_csv(BASE / "results/covalent_strategy_comparison.csv")
hpa = pd.read_csv(BASE / "results/hpa_tissue_expression.csv")
mut = pd.read_csv(BASE / "results/tcga_mutation_frequencies.csv")
cons = pd.read_csv(BASE / "results/cross_species_conservation.csv")
ot = pd.read_csv(BASE / "results/opentargets_disease.csv")

# ============================================================
# Figure 6: Cancer-Type Specific Expression Heatmap
# ============================================================
print("Generating Figure 6...")

# Prepare heatmap data: top 15 genes × all cancer types
top15_genes = enhanced.head(15)["gene"].tolist()
top15_genes = [g for g in top15_genes if g in ct_df["gene"].values]
ct_subset = ct_df[ct_df["gene"].isin(top15_genes)].copy()

# Build matrix
cancers = sorted(set(ct_subset["best_cancer"].tolist() + ct_subset["second_cancer"].tolist()))
# Also get all unique cancers from the original data
all_cancers = ["Colorectal","Lung Adeno","Pancreatic","Prostate","Thyroid",
               "Breast","Bladder","Kidney Clear Cell","Kidney Papillary",
               "Ovarian","Stomach","Liver","Head&Neck","Melanoma","Uterine",
               "Esophageal","Cervical","Sarcoma","Glioblastoma","Lung Squamous"]

# For each gene, we have best_cancer and its log2tpm, plus second
heatmap_data = pd.DataFrame(index=top15_genes, columns=all_cancers)

for _, row in ct_subset.iterrows():
    gene = row["gene"]
    if gene in heatmap_data.index:
        # Best cancer
        if row["best_cancer"] in heatmap_data.columns:
            heatmap_data.loc[gene, row["best_cancer"]] = row["best_log2tpm"]
        # Second cancer
        if row["second_cancer"] in heatmap_data.columns:
            heatmap_data.loc[gene, row["second_cancer"]] = row["second_log2tpm"]

heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')
# Keep only cancers with data
heatmap_data = heatmap_data.dropna(axis=1, how='all')

fig6, ax = plt.subplots(figsize=(14, 6))
sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd',
            linewidths=0.5, ax=ax, cbar_kws={'label': 'Median log2(TPM)'},
            vmin=0, vmax=12)
ax.set_title('Cancer-Type Specific Expression of Top LDT Targets (log2 TPM)')
ax.set_xlabel('Cancer Type')
ax.set_ylabel('Gene')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
fig6.savefig(FIG_DIR / "figure6_cancer_heatmap.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 6 saved")

# ============================================================
# Figure 7: Covalent Strategy Comparison
# ============================================================
print("Generating Figure 7...")

# Count targets best suited for each strategy
strategies = ["LDT_NAS", "SuFEx_CTR", "Traditional_Acrylamide"]
if "best_strategy" in covalent.columns:
    strat_counts = covalent["best_strategy"].value_counts()

fig7, axes = plt.subplots(1, 2, figsize=(14, 5))

# (a) Strategy preference pie chart
ax = axes[0]
labels = {"LDT_NAS": "LDT-NASA\n(Lys primary)",
          "SuFEx_CTR": "SuFEx-CTR\n(Tyr primary)",
          "Traditional_Acrylamide": "Traditional\n(Cys primary)"}
colors_strat = {"LDT_NAS": "#e74c3c", "SuFEx_CTR": "#3498db", "Traditional_Acrylamide": "#2ecc71"}
valid_counts = {k: v for k, v in strat_counts.items() if k in labels}
ax.pie(valid_counts.values(), labels=[labels[k] for k in valid_counts.keys()],
       autopct='%1.1f%%', colors=[colors_strat[k] for k in valid_counts.keys()],
       startangle=90)
ax.set_title('Optimal Covalent Strategy per Target')

# (b) Score distribution comparison
ax = axes[1]
gene_examples = enhanced.head(15)["gene"].tolist()
for i, sk in enumerate(strategies):
    col = f"{sk}_score"
    if col in covalent.columns:
        vals = covalent[col].dropna()
        ax.hist(vals, bins=20, alpha=0.5, label=labels.get(sk, sk), color=colors_strat.get(sk, '#95a5a6'))
ax.set_xlabel('Covalent Suitability Score')
ax.set_ylabel('Number of Targets')
ax.set_title('Score Distribution: Three Covalent Strategies')
ax.legend(frameon=False)

plt.tight_layout()
fig7.savefig(FIG_DIR / "figure7_covalent_comparison.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 7 saved")

# ============================================================
# Figure 8: Tool Comparison Radar + Multi-Dimension Validation
# ============================================================
print("Generating Figure 8...")

fig8, axes = plt.subplots(1, 2, figsize=(14, 5))

# (a) Bar chart comparing feature counts across tools
ax = axes[0]
tools = ["LDT-TargetDB", "DrugMap", "ImmunoTar", "TCSA"]
features_count = [19, 11, 9, 8]
colors_bar = ['#e74c3c', '#3498db', '#2ecc71', '#95a5a6']
bars = ax.bar(tools, features_count, color=colors_bar, edgecolor='white', width=0.6)
for bar, count in zip(bars, features_count):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            str(count), ha='center', fontweight='bold', fontsize=12)
ax.set_ylabel('Number of Features (of 20)')
ax.set_title('Feature Completeness Comparison')
ax.set_ylim(0, 22)

# (b) Validation summary: integrate HPA, mutations, conservation, disease
ax = axes[1]
# Radiopharmaceutical suitability radar
# For top 6 targets, show their multi-dimensional scores
top6 = enhanced.head(6)
genes_radar = top6["gene"].tolist()

# Build a summary table as bar chart
categories = ['TSI High\n(>0.7)', 'LDT Score\n(>0.5)', 'Low Mutation\n(<30)', 'Conserved\n(>5 models)', 'Cancer-Assoc\n(OpenTargets)']
n_qualifying = [
    (enhanced["tsi_norm"] > 0.7).sum(),
    (enhanced["ldt_norm"] > 0.5).sum(),
    (mut["total_mutations"] < 30).sum() if "total_mutations" in mut.columns else 0,
    (cons["conservation_level"] == "High (>5)").sum() if "conservation_level" in cons.columns else 0,
    (ot["n_cancer_associations"] > 0).sum() if "n_cancer_associations" in ot.columns else 0,
]
ax.barh(categories, n_qualifying, color='#e74c3c', height=0.6, alpha=0.8)
for i, (cat, val) in enumerate(zip(categories, n_qualifying)):
    ax.text(val + 1, i, f"{val}", va='center', fontweight='bold')
ax.set_xlabel('Number of Qualifying Targets (of 195)')
ax.set_title('Multi-Dimensional Validation Summary')
ax.set_xlim(0, max(n_qualifying) * 1.3)

plt.tight_layout()
fig8.savefig(FIG_DIR / "figure8_tool_comparison.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 8 saved")

# ============================================================
# Tables
# ============================================================
print("Generating Tables...")

# Table 1: Database Statistics
table1 = pd.DataFrame({
    "Category": [
        "Surface proteins (SURFY)", "With TCGA expression data", "With GTEx normal data",
        "AlphaFold structures downloaded", "With detectable pockets", "With LDT-suitable nucleophiles",
        "Cancer types analyzed", "TCGA tumor samples", "GTEx normal samples",
        "Cell lines (DepMap)", "GO terms enriched (TOP30)", "STRING PPI edges (TOP30)",
        "HPA protein IHC entries queried", "Open Targets disease queries", "TCGA cancer types (mutation)",
    ],
    "Count": [
        "2,799", "2,663", "2,663",
        "2,490", "~2,480 (99.6%)", "~2,360 (95%)",
        "32", "9,186", "7,862",
        "1,208", "20 (p<0.05)", "45",
        "1,199,675", "30 genes", "37",
    ]
})
table1.to_csv(BASE / "results" / "table1_database_statistics.csv", index=False)
print("  Table 1 saved")

# Table 3: TOP20 Target Details
cols_top20 = ["gene", "enhanced_score", "tsi_norm", "ldt_norm", "top_res_type",
              "n_nucleophiles", "depmap_label", "mean_plddt", "n_cancers_expressed"]
cols_avail = [c for c in cols_top20 if c in enhanced.columns]

# Merge with cancer-type and additional data
table3 = enhanced.head(20)[cols_avail].copy()
if "gene" in ct_df.columns:
    table3 = table3.merge(ct_df[["gene","best_cancer","best_log2tpm","n_cancers_positive"]],
                          on="gene", how="left")
table3.to_csv(BASE / "results" / "table3_top20_targets.csv", index=False)
print("  Table 3 saved")

# Table 4: Known Targets Validation
table4 = known[["gene", "standard_name", "tsi", "tsi_rank", "rank", "final_score",
                "n_pockets", "ldt_score", "top_res_type"]].copy()
table4.columns = ["Gene", "Standard Name", "TSI", "TSI Rank", "LDT Rank",
                  "LDT Score", "Pockets", "Raw LDT", "Best Residue"]
table4.to_csv(BASE / "results" / "table4_known_target_validation.csv", index=False)
print("  Table 4 saved")

# Table 2: Tool Comparison (already exists, just copy)
print("  Table 2 (tool_comparison.csv) already exists")

print(f"\nAll figures and tables saved to: {FIG_DIR} and {BASE / 'results'}")
print("Done.")
