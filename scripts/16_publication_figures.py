"""
Publication-ready figures for NAR submission
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.lines import Line2D
import seaborn as sns

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
FIG_DIR = BASE / "results/figures"
FIG_DIR.mkdir(exist_ok=True)

# ===== PROFESSIONAL STYLE =====
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': False,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
})

# Consistent color palette
C_BLUE = '#2471A3'
C_RED = '#CB4335'
C_GREEN = '#27AE60'
C_ORANGE = '#E67E22'
C_PURPLE = '#8E44AD'
C_TEAL = '#17A589'
C_GRAY = '#7F8C8D'
C_DARK = '#2C3E50'

# Residue colors (colorblind-friendly)
RES_COLORS = {'LYS': '#E74C3C', 'CYS': '#F39C12', 'TYR': '#2980B9', 'SER': '#95A5A6', 'HIS': '#27AE60'}

def panel_label(ax, text, x=-0.06, y=1.02):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# Load data
enhanced = pd.read_csv(BASE / "results/enhanced_final_ranking.csv")
tsi = pd.read_csv(BASE / "data/expression/tsi_results.csv")
full = pd.read_csv(BASE / "results/full_proteome_ranking.csv")
ct = pd.read_csv(BASE / "results/cancer_type_specific_expression.csv")
covalent = pd.read_csv(BASE / "results/covalent_strategy_comparison.csv")

top30 = enhanced.head(30)
top20 = enhanced.head(20)

# ============================================================
# Figure 1: TSI Overview + Pocket + TSI vs LDT
# ============================================================
print("Figure 1...")
fig1, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# (A) TSI histogram
ax = axes[0]
ax.hist(tsi["tsi"], bins=90, color=C_BLUE, alpha=0.75, edgecolor='white', linewidth=0.3)
ax.axvline(x=0, color=C_RED, linestyle='--', linewidth=1.2)
ax.text(0.1, ax.get_ylim()[1]*0.9, 'TSI = 0', color=C_RED, fontsize=8)
ax.set_xlabel('Tumor Specificity Index (TSI)')
ax.set_ylabel('Number of Proteins')
panel_label(ax, 'A')

# (B) Pocket histogram
ax = axes[1]
counts = full["n_pockets"].value_counts().sort_index()
ax.bar(counts.index[:40], counts.values[:40], color=C_GREEN, alpha=0.8, edgecolor='white', linewidth=0.3)
ax.axvline(x=full["n_pockets"].median(), color=C_RED, linestyle='--', linewidth=1.2,
           label=f'Median: {full["n_pockets"].median():.0f}')
ax.legend(frameon=False, fontsize=7)
ax.set_xlabel('Number of Pockets')
ax.set_ylabel('Count')
panel_label(ax, 'B')

# (C) TSI vs LDT scatter
ax = axes[2]
for rt in ['LYS','CYS','TYR','SER']:
    s = enhanced[enhanced["top_res_type"]==rt]
    ax.scatter(s["tsi_norm"], s["ldt_norm"], c=RES_COLORS[rt], label=rt, alpha=0.5, s=12, edgecolors='none')
ax.set_xlabel('TSI (normalized)')
ax.set_ylabel('LDT Score (normalized)')
ax.legend(frameon=False, title='Best Residue', title_fontsize=8, fontsize=7, ncol=2, loc='lower right')
panel_label(ax, 'C')

plt.tight_layout(pad=1.5)
fig1.savefig(FIG_DIR / "figure1_overview.png")
plt.close()

# ============================================================
# Figure 2: TOP20 Ranking + Nucleophile Pie
# ============================================================
print("Figure 2...")
fig2, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# (A) Horizontal bar chart
ax = axes[0]
genes = top20["gene"].tolist()[::-1]
scores = top20["final_score"].tolist()[::-1]
types = top20["top_res_type"].tolist()[::-1]
bars = ax.barh(range(len(genes)), scores, height=0.7)
for i, (g, s, t) in enumerate(zip(genes, scores, types)):
    bars[i].set_color(RES_COLORS.get(t, C_GRAY))
    bars[i].set_alpha(0.85)
    ax.text(s+0.01, i, g, fontsize=7.5, va='center', fontweight='bold')
ax.set_xlim(0, max(scores)*1.35)
ax.set_yticks([])
ax.set_xlabel('Composite LDT Score')
panel_label(ax, 'A')

# Legend
patches = [Patch(color=c, label=r) for r,c in RES_COLORS.items() if r in types]
ax.legend(handles=patches, frameon=False, fontsize=7, ncol=4, loc='lower right')

# (B) Pie chart
ax = axes[1]
rc = full["top_res_type"].value_counts()
colors_pie = [RES_COLORS.get(r, C_GRAY) for r in rc.index]
wedges, texts, autotexts = ax.pie(rc.values, labels=rc.index, autopct='%1.1f%%',
    colors=colors_pie, startangle=90, pctdistance=0.75,
    wedgeprops=dict(width=0.4, edgecolor='white', linewidth=0.5))
for t in autotexts: t.set_fontsize(8)
for t in texts: t.set_fontsize(9)
ax.set_title('Nucleophile Distribution (n={})'.format(len(full)), fontsize=10, pad=12)
panel_label(ax, 'B')

plt.tight_layout(pad=1.5)
fig2.savefig(FIG_DIR / "figure2_top20.png")
plt.close()

# ============================================================
# Figure 3: Known Targets Validation + DepMap
# ============================================================
print("Figure 3...")
fig3, axes = plt.subplots(1, 2, figsize=(13, 4.5))

# (A) TSI histogram + known target lines
ax = axes[0]
ax.hist(tsi["tsi"], bins=100, color=C_BLUE, alpha=0.6, edgecolor='white', linewidth=0.2)

known = ["TACSTD2","MET","LY6E","FOLH1","ERBB2","FAP","DLL3","CLDN6","SSTR2","EGFR"]
colors_known = ["#E74C3C","#27AE60","#2980B9","#1ABC9C","#F39C12","#8E44AD","#E67E22","#16A085","#C0392B","#D35400"]
right_edge = ax.get_xlim()[1] * 0.82

for i, (gene, color) in enumerate(zip(known, colors_known)):
    r = tsi[tsi["gene"].str.upper()==gene.upper()]
    if len(r)==0: continue
    tv = r["tsi"].values[0]
    ax.axvline(x=tv, color=color, linewidth=1.8, linestyle='--', alpha=0.8)
    y_pos = 0.94 - i*0.09
    ax.text(right_edge, y_pos, gene, fontsize=7, color=color, fontweight='bold',
            ha='left', va='center', transform=ax.get_xaxis_transform())
    ax.plot([tv, right_edge], [y_pos, y_pos], color=color, alpha=0.2, linewidth=0.5,
            transform=ax.get_xaxis_transform())

ax.set_xlabel('Tumor Specificity Index (TSI)')
ax.set_ylabel('Count')
panel_label(ax, 'A')

# (B) DepMap scatter
ax = axes[1]
dm = enhanced[enhanced["depmap_score"].notna()]
sc = ax.scatter(dm["tsi_norm"], dm["depmap_score"], c=dm["ldt_norm"],
                cmap='RdYlGn', alpha=0.5, s=15, edgecolors='none')
ax.axhline(y=-0.5, color=C_RED, linestyle='--', linewidth=1, alpha=0.6)
ax.text(1.0, -0.46, 'Essential', ha='right', fontsize=7, color=C_RED, transform=ax.get_xaxis_transform())
cbar = plt.colorbar(sc, ax=ax, aspect=30, pad=0.02)
cbar.set_label('LDT Score', fontsize=8)
cbar.ax.tick_params(labelsize=7)
ax.set_xlabel('TSI (normalized)')
ax.set_ylabel('DepMap Chronos Score')
panel_label(ax, 'B')

plt.tight_layout(pad=1.5)
fig3.savefig(FIG_DIR / "figure3_validation.png")
plt.close()

# ============================================================
# Figure 4: Expression Features + Cancer Heatmap
# ============================================================
print("Figure 4...")
fig4, axes = plt.subplots(1, 2, figsize=(13, 5))

# (A) Expression levels
ax = axes[0]
t30 = enhanced.head(30)[::-1]
ec = "tumor_median_log2" if "tumor_median_log2" in t30.columns else "tsi_norm"
ax.barh(range(30), t30[ec], height=0.65, color=C_BLUE, alpha=0.8)
ax.set_yticks(range(30))
ax.set_yticklabels(t30["gene"].values, fontsize=7)
ax.set_xlabel('Tumor Median log2(TPM)')
panel_label(ax, 'A')

# (B) Cancer heatmap (top 15 genes x cancer types)
ax = axes[1]
top15 = enhanced.head(15)["gene"].tolist()
top15_in_ct = [g for g in top15 if g in ct["gene"].values]
ct_sub = ct[ct["gene"].isin(top15_in_ct)]

all_cancers = ["Colorectal","Lung Adeno","Pancreatic","Prostate","Thyroid",
               "Breast","Bladder","Kidney Clear Cell","Kidney Papillary",
               "Ovarian","Stomach","Liver","Head&Neck","Melanoma","Uterine",
               "Esophageal","Cervical","Sarcoma","Glioblastoma","Lung Squamous"]
hdata = pd.DataFrame(index=top15_in_ct, columns=all_cancers)
for _, row in ct_sub.iterrows():
    g = row["gene"]
    if g in hdata.index:
        if row["best_cancer"] in hdata.columns:
            hdata.loc[g, row["best_cancer"]] = row["best_log2tpm"]
        if row["second_cancer"] in hdata.columns:
            hdata.loc[g, row["second_cancer"]] = row["second_log2tpm"]
hdata = hdata.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')

sns.heatmap(hdata, annot=True, fmt='.1f', cmap='YlOrRd', linewidths=0.3,
            ax=ax, cbar_kws={'label': 'log2 TPM', 'shrink': 0.7}, annot_kws={'fontsize':6})
ax.set_xlabel('Cancer Type', fontsize=8)
ax.set_ylabel('Gene', fontsize=8)
ax.tick_params(labelsize=7)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
panel_label(ax, 'B')

plt.tight_layout(pad=1.5)
fig4.savefig(FIG_DIR / "figure4_cancer_types.png")
plt.close()

# ============================================================
# Figure 5: pLDDT Quality + Pocket Count
# ============================================================
print("Figure 5...")
fig5, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# (A) pLDDT vs LDT
ax = axes[0]
valid = enhanced.dropna(subset=["mean_plddt","ldt_norm"])
for rt in ['LYS','CYS','TYR','SER']:
    s = valid[valid["top_res_type"]==rt]
    ax.scatter(s["mean_plddt"], s["ldt_norm"], c=RES_COLORS[rt], label=rt, alpha=0.4, s=10, edgecolors='none')

# Top 10 labels with arrows
top10 = enhanced.head(10).dropna(subset=["mean_plddt"])
lx = ax.get_xlim()[0] + 1.5
ly = ax.get_ylim()[1] * 0.95
sp = (ax.get_ylim()[1]-ax.get_ylim()[0])*0.045
text_objs = []
for i, (_, row) in enumerate(top10.iterrows()):
    t = ax.text(lx, ly-i*sp, row["gene"], fontsize=6.5, color=RES_COLORS.get(row["top_res_type"],'#333'),
                fontweight='bold', ha='left', va='center')
    text_objs.append((t, row))

fig5.canvas.draw()
for t, row in text_objs:
    bbox = t.get_window_extent(renderer=fig5.canvas.get_renderer())
    bbox_d = ax.transData.inverted().transform(bbox)
    ax.annotate('', xy=(row["mean_plddt"], row["ldt_norm"]),
                xytext=(bbox_d[1,0], (bbox_d[1,1]+bbox_d[0,1])/2),
                arrowprops=dict(arrowstyle='->', color=RES_COLORS.get(row["top_res_type"],'#333'), lw=0.6, alpha=0.4))

ax.set_xlabel('Mean pLDDT')
ax.set_ylabel('LDT Score (normalized)')
ax.legend(frameon=False, fontsize=7, ncol=2)
panel_label(ax, 'A')

# (B) Pocket distribution
ax = axes[1]
ax.hist(full["n_pockets"], bins=60, color=C_PURPLE, alpha=0.8, edgecolor='white', linewidth=0.3)
ax.axvline(x=full["n_pockets"].median(), color=C_RED, linestyle='--', linewidth=1)
ax.text(full["n_pockets"].median()+1, ax.get_ylim()[1]*0.9,
        f'Median: {full["n_pockets"].median():.0f}', fontsize=8, color=C_RED)
ax.set_xlabel('Number of Pockets')
ax.set_ylabel('Count')
panel_label(ax, 'B')

plt.tight_layout(pad=1.5)
fig5.savefig(FIG_DIR / "figure5_quality.png")
plt.close()

# ============================================================
# Figure 6: Covalent Strategy Comparison
# ============================================================
print("Figure 6...")
fig6, axes = plt.subplots(1, 2, figsize=(12, 4.5))

strat_labels = {"LDT_NAS": "LDT-NASA\n(Lys primary)",
                "SuFEx_CTR": "SuFEx-CTR\n(Tyr primary)",
                "Traditional_Acrylamide": "Traditional\n(Cys primary)"}
strat_colors = {"LDT_NAS": "#E74C3C", "SuFEx_CTR": "#2980B9", "Traditional_Acrylamide": "#27AE60"}

# (A) Pie
ax = axes[0]
if "best_strategy" in covalent.columns:
    sc = covalent["best_strategy"].value_counts()
    vc = {k: v for k,v in sc.items() if k in strat_labels}
    ax.pie(vc.values(), labels=[strat_labels[k] for k in vc.keys()],
           autopct='%1.1f%%', colors=[strat_colors[k] for k in vc.keys()],
           startangle=90, textprops={'fontsize':8})
panel_label(ax, 'A')

# (B) Score distribution
ax = axes[1]
for sk in strat_labels:
    col = f"{sk}_score"
    if col in covalent.columns:
        ax.hist(covalent[col].dropna(), bins=30, alpha=0.5, label=strat_labels[sk],
                color=strat_colors[sk], edgecolor='white', linewidth=0.3)
ax.set_xlabel('Covalent Suitability Score')
ax.set_ylabel('Number of Targets')
ax.legend(frameon=False, fontsize=7)
panel_label(ax, 'B')

plt.tight_layout(pad=1.5)
fig6.savefig(FIG_DIR / "figure6_covalent.png")
plt.close()

# ============================================================
# Figure 7: Feature Comparison + Validation Summary
# ============================================================
print("Figure 7...")
fig7, axes = plt.subplots(1, 2, figsize=(13, 4.5))

# (A) Tool comparison bar
ax = axes[0]
tools = ["LDT-TargetDB","DrugMap","ImmunoTar","TCSA"]
features = [19,11,9,8]
bar_colors = [C_RED, C_BLUE, C_GREEN, C_GRAY]
bars = ax.bar(tools, features, color=bar_colors, width=0.55, edgecolor='white', linewidth=0.5)
for b, c in zip(bars, features):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.3, str(c), ha='center', fontweight='bold', fontsize=13, color='#333')
ax.set_ylabel('Features (of 20)')
ax.set_ylim(0, 23)
ax.tick_params(axis='x', labelsize=9)
panel_label(ax, 'A')

# (B) Validation summary
ax = axes[1]
hpa_data = pd.read_csv(BASE / "results/hpa_tissue_expression.csv")
mut_data = pd.read_csv(BASE / "results/tcga_mutation_frequencies.csv")
cons_data = pd.read_csv(BASE / "results/cross_species_conservation.csv")
ot_data = pd.read_csv(BASE / "results/opentargets_disease.csv")

cats = ['TSI High\n(>0.7)','LDT Score\n(>0.5)','Low Mutation\n(<30)','Conserved\n(>5 models)','Cancer-Assoc\n(OpenTargets)']
vals = [
    (enhanced["tsi_norm"]>0.7).sum(),
    (enhanced["ldt_norm"]>0.5).sum(),
    (mut_data["total_mutations"]<30).sum() if "total_mutations" in mut_data.columns else 0,
    (cons_data["conservation_level"]=="High (>5)").sum() if "conservation_level" in cons_data.columns else 0,
    (ot_data["n_cancer_associations"]>0).sum() if "n_cancer_associations" in ot_data.columns else 0,
]
ax.barh(cats, vals, color=C_RED, height=0.55, alpha=0.85)
for i, (c, v) in enumerate(zip(cats, vals)):
    ax.text(v+2, i, str(v), va='center', fontweight='bold', fontsize=10)
ax.set_xlabel('Number of Qualifying Targets')
ax.set_xlim(0, max(vals)*1.35)
panel_label(ax, 'B')

plt.tight_layout(pad=1.5)
fig7.savefig(FIG_DIR / "figure7_tools.png")
plt.close()

print("All 7 figures saved at 600 DPI")
