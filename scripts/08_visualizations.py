"""
论文可视化：Figure 1-5生成脚本
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
FIG_DIR = BASE / "results" / "figures"
FIG_DIR.mkdir(exist_ok=True)

# 全局样式
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
})

# ============================================================
# Figure 1: 分析流程概况 + TSI分布
# ============================================================
print("Generating Figure 1...")
enhanced = pd.read_csv(BASE / "results" / "enhanced_final_ranking.csv")
tsi = pd.read_csv(BASE / "data" / "expression" / "tsi_results.csv")

fig1, axes = plt.subplots(1, 3, figsize=(14, 4))

# (a) TSI distribution
ax = axes[0]
tsi_valid = tsi[tsi["tsi"].notna() & (tsi["tsi"] > -5)]
ax.hist(tsi_valid["tsi"], bins=80, color='#3498db', edgecolor='white', alpha=0.8)
ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='TSI=0')
ax.set_xlabel('Tumor Specificity Index (TSI)')
ax.set_ylabel('Number of Surface Proteins')
ax.set_title('TSI Distribution (2,663 surface proteins)')
ax.legend(frameon=False)

# (b) Pocket distribution
ax = axes[1]
pocket_counts = enhanced["n_pockets"].value_counts().sort_index()
ax.bar(pocket_counts.index[:50], pocket_counts.values[:50], color='#2ecc71', edgecolor='white')
ax.set_xlabel('Number of Pockets per Protein')
ax.set_ylabel('Number of Proteins')
ax.set_title('Pocket Detection (n=195)')

# (c) LDT score vs TSI scatter
ax = axes[2]
colors = {"LYS": "#e74c3c", "CYS": "#f39c12", "TYR": "#3498db", "SER": "#95a5a6", "HIS": "#2ecc71"}
for rt in ["LYS", "CYS", "TYR", "SER"]:
    subset = enhanced[enhanced["top_res_type"] == rt]
    if len(subset) > 0:
        ax.scatter(subset["tsi_norm"], subset["ldt_norm"], c=colors[rt], label=rt,
                   alpha=0.7, s=40, edgecolors='white')
ax.set_xlabel('TSI (normalized)')
ax.set_ylabel('LDT Transferability Score')
ax.set_title('TSI vs LDT by Nucleophile Type')
ax.legend(frameon=False, title='Best Residue')

plt.tight_layout()
fig1.savefig(FIG_DIR / "figure1_overview.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 1 saved")

# ============================================================
# Figure 2: TOP20靶点瀑布图 + 残基类型分布
# ============================================================
print("Generating Figure 2...")
top20 = enhanced.head(20).copy()
top20 = top20.iloc[::-1]  # Reverse for horizontal bar

fig2, axes = plt.subplots(1, 2, figsize=(14, 6))

# (a) Waterfall/bar of composite score
ax = axes[0]
ax.barh(range(len(top20)), top20["final_score"], color='#2c3e50', height=0.7)
# Color-code by residue
for i, (_, row) in enumerate(top20.iterrows()):
    c = colors.get(row["top_res_type"], "#2c3e50")
    ax.barh(i, row["final_score"], color=c, height=0.7, alpha=0.8)
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20["gene"].values)
ax.set_xlabel('Composite LDT Score')
ax.set_title('Top 20 LDT Radiopharmaceutical Targets')

# Legend
patches = [mpatches.Patch(color=c, label=rt) for rt, c in colors.items() if any(top20["top_res_type"] == rt)]
ax.legend(handles=patches, frameon=False, title='Nucleophile')

# (b) Nucleophile type distribution (all 195)
ax = axes[1]
res_counts = enhanced["top_res_type"].value_counts()
explode = [0.05 if i == 0 else 0 for i in range(len(res_counts))]
wedges, texts, autotexts = ax.pie(
    res_counts.values, labels=res_counts.index, autopct='%1.1f%%',
    colors=[colors.get(r, '#95a5a6') for r in res_counts.index],
    explode=explode, startangle=90
)
ax.set_title('Nucleophile Distribution in Best Pockets')

plt.tight_layout()
fig2.savefig(FIG_DIR / "figure2_top20.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 2 saved")

# ============================================================
# Figure 3: 已知靶点验证 + DepMap
# ============================================================
print("Generating Figure 3...")
known = pd.read_csv(BASE / "results" / "known_targets_validation.csv")

# Only known targets with data
# For now, use TSI data
known_tsi = tsi[tsi["gene"].str.upper().isin(
    [g.upper() for g in ["FOLH1","FAP","SSTR2","TACSTD2","MET","LY6E","ERBB2","EGFR","DLL3","CLDN6"]]
)].copy()
known_tsi["label"] = known_tsi["gene"].str.upper()

fig3, axes = plt.subplots(1, 2, figsize=(14, 5))

# (a) Known targets in our ranking
ax = axes[0]
# All genes TSI
all_tsi_sorted = tsi.sort_values("tsi", ascending=False).reset_index(drop=True)
all_tsi_sorted["rank_pct"] = range(1, len(all_tsi_sorted)+1)

for _, row in known_tsi.head(10).iterrows():
    rank_pos = tsi.sort_values("tsi", ascending=False).reset_index(drop=True)
    gene_str = str(row["gene"]).upper()
    r = rank_pos[rank_pos["gene"].str.upper() == gene_str].index
    if len(r) > 0:
        rank_pct = r[0] / len(all_tsi_sorted) * 100
        res_type = '#bdc3c7'
        matched = enhanced[enhanced["gene"].str.upper() == str(row["gene"]).upper()]
        if len(matched) > 0:
            res_type = colors.get(matched["top_res_type"].values[0], '#bdc3c7')
        ax.axvline(x=rank_pct, color=res_type, alpha=0.5, linewidth=1)
        ax.text(rank_pct+0.5, 0.5, row["gene"], rotation=90, fontsize=7, va='center')

ax.hist(np.clip(tsi["tsi"], -2, 3), bins=100, color='#3498db', alpha=0.6, edgecolor='white')
ax.set_xlabel('Tumor Specificity Index (TSI)')
ax.set_ylabel('Count')
ax.set_title('Known Nuclear Medicine Targets in TSI Landscape')

# (b) DepMap comparison
ax = axes[1]
depmap_data = enhanced[enhanced["depmap_score"].notna()].copy()
ax.scatter(depmap_data["tsi_norm"], depmap_data["depmap_score"],
           c=depmap_data["ldt_norm"], cmap='RdYlGn', alpha=0.6, s=30)
ax.axhline(y=-0.5, color='red', linestyle='--', alpha=0.5, label='Essential threshold')
ax.set_xlabel('TSI (normalized)')
ax.set_ylabel('DepMap Chronos Score')
ax.set_title('Gene Essentiality vs Tumor Specificity')
ax.legend(frameon=False)

plt.tight_layout()
fig3.savefig(FIG_DIR / "figure3_validation.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 3 saved")

# ============================================================
# Figure 4: 癌种特异性热图
# ============================================================
print("Generating Figure 4...")
# 使用肿瘤表达特征（替代癌种特异性，因为TCGA barcode不含癌种缩写）
fig4, axes = plt.subplots(1, 2, figsize=(14, 5))

# (a) Tumor expression levels of top 30 targets
ax = axes[0]
top30 = enhanced.head(30).copy()
expr_col = "tumor_median_log2" if "tumor_median_log2" in top30.columns else "tsi_norm"
expr_vals = top30[expr_col] if expr_col in top30.columns else top30["tsi_norm"]
ax.barh(range(len(top30)), expr_vals,
        color='#3498db', height=0.7, alpha=0.8)
ax.set_yticks(range(len(top30)))
ax.set_yticklabels(top30["gene"].values, fontsize=8)
ax.set_xlabel('Tumor Median log2(TPM)')
ax.set_title('Tumor Expression Levels of Top 30 Targets')
ax.invert_yaxis()

# (b) Tumor positive rate distribution
ax = axes[1]
pos_col = "tumor_positive_rate" if "tumor_positive_rate" in top30.columns else "tsi_norm"
pos_vals = top30[pos_col] if pos_col in top30.columns else top30["tsi_norm"]
ax.barh(range(len(top30)), pos_vals,
        color='#e74c3c', height=0.7, alpha=0.8)
ax.set_yticks(range(len(top30)))
ax.set_yticklabels([])
ax.set_xlabel('Tumor Positive Rate')
ax.set_title('Expression Coverage Across Tumors')
ax.invert_yaxis()

plt.tight_layout()
fig4.savefig(FIG_DIR / "figure4_expression_features.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 4 saved")

# ============================================================
# Figure 5: pLDDT置信度 + 口袋可药性
# ============================================================
print("Generating Figure 5...")
plddt = pd.read_csv(BASE / "results" / "plddt_confidence.csv")

fig5, axes = plt.subplots(1, 2, figsize=(12, 4))

# (a) pLDDT vs LDT score
ax = axes[0]
valid = enhanced.dropna(subset=["mean_plddt", "ldt_norm"])
sc = ax.scatter(valid["mean_plddt"], valid["ldt_norm"],
                c=valid["tsi_norm"], cmap='viridis', alpha=0.7, s=40, edgecolors='white')
ax.set_xlabel('Mean pLDDT')
ax.set_ylabel('LDT Score (normalized)')
ax.set_title('Structure Confidence vs LDT Transferability')
plt.colorbar(sc, ax=ax, label='TSI')

# Highlight top targets
for _, row in enhanced.head(10).dropna(subset=["mean_plddt"]).iterrows():
    ax.annotate(row["gene"], (row["mean_plddt"], row["ldt_norm"]),
                fontsize=7, alpha=0.8, xytext=(5, 5), textcoords='offset points')

# (b) Pocket volume distribution (from full ranking if available)
ax = axes[1]
vol_col = "max_pocket_volume" if "max_pocket_volume" in enhanced.columns else "n_pockets"
vol_data = enhanced[vol_col].dropna() if vol_col in enhanced.columns else enhanced["n_pockets"].dropna()
ax.hist(vol_data, bins=min(50, len(vol_data.unique())), color='#9b59b6', edgecolor='white', alpha=0.8)
ax.axvline(x=vol_data.median(), color='red', linestyle='--',
           label=f'Median: {vol_data.median():.1f}')
ax.set_xlabel(vol_col.replace('_', ' ').title())
ax.set_ylabel('Count')
ax.set_title('Pocket Metric Distribution')
ax.legend(frameon=False)

plt.tight_layout()
fig5.savefig(FIG_DIR / "figure5_quality.png", dpi=200, bbox_inches='tight')
plt.close()
print("  Figure 5 saved")

print(f"\nAll figures saved to: {FIG_DIR}")
print("Done.")
