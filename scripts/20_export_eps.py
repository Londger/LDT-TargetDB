import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np, seaborn as sns
from pathlib import Path
from PIL import Image
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
FIG = BASE / "results/figures/eps"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9,
    'axes.titlesize': 10, 'savefig.bbox': 'tight',
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})

RC = {'LYS': '#E74C3C', 'CYS': '#F39C12', 'TYR': '#2980B9', 'SER': '#95A5A6'}
BL, RE, GR, PU = '#2471A3', '#CB4335', '#27AE60', '#8E44AD'

e = pd.read_csv(BASE / 'results/enhanced_final_ranking.csv')
t = pd.read_csv(BASE / 'data/expression/tsi_results.csv')
f = pd.read_csv(BASE / 'results/full_proteome_ranking.csv')
ct = pd.read_csv(BASE / 'results/cancer_type_specific_expression.csv')

def save(fig, name):
    fig.savefig(str(FIG / name), format='eps')
    print('  ' + name + ' saved')

# Fig 1
print('Fig1...')
fig1, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax[0].hist(t['tsi'], bins=90, color=BL, alpha=0.75, edgecolor='white', linewidth=0.3)
ax[0].axvline(x=0, color=RE, linestyle='--', linewidth=1.2)
ax[0].set_xlabel('TSI')
ax[0].set_ylabel('Count')
counts = f['n_pockets'].value_counts().sort_index()
ax[1].bar(counts.index[:40], counts.values[:40], color=GR, alpha=0.8, edgecolor='white', linewidth=0.3)
ax[1].set_xlabel('Pockets')
ax[1].set_ylabel('Count')
for rt in ['LYS', 'CYS', 'TYR', 'SER']:
    s = e[e['top_res_type'] == rt]
    ax[2].scatter(s['tsi_norm'], s['ldt_norm'], c=RC[rt], label=rt, alpha=0.5, s=12, edgecolors='none')
ax[2].set_xlabel('TSI (norm)')
ax[2].set_ylabel('LDT (norm)')
ax[2].legend(fontsize=7, ncol=2, loc='lower right')
plt.tight_layout(pad=1.5)
save(fig1, 'figure1_overview.eps')
plt.close()

# Fig 2
print('Fig2...')
fig2, ax = plt.subplots(1, 2, figsize=(13, 5.5))
t20 = e.head(20)
genes = t20['gene'].tolist()[::-1]
scores = t20['final_score'].tolist()[::-1]
types = t20['top_res_type'].tolist()[::-1]
bars = ax[0].barh(range(len(genes)), scores, height=0.7)
for i in range(len(genes)):
    bars[i].set_color(RC.get(types[i], '#7F8C8D'))
    bars[i].set_alpha(0.85)
    ax[0].text(scores[i] + 0.01, i, genes[i], fontsize=7.5, va='center', fontweight='bold')
ax[0].set_xlim(0, max(scores) * 1.35)
ax[0].set_yticks([])
ax[0].set_xlabel('Composite Score')
rc = f['top_res_type'].value_counts()
ax[1].pie(rc.values, labels=rc.index, autopct='%1.1f%%',
    colors=[RC.get(r, '#7F8C8D') for r in rc.index],
    startangle=90, pctdistance=0.75,
    wedgeprops=dict(width=0.4, edgecolor='white', linewidth=0.5))
plt.tight_layout(pad=1.5)
save(fig2, 'figure2_top20.eps')
plt.close()

# Fig 3
print('Fig3...')
fig3, ax = plt.subplots(1, 2, figsize=(13, 4.5))
ax[0].hist(t['tsi'], bins=100, color=BL, alpha=0.6, edgecolor='white', linewidth=0.2)
kn = ['TACSTD2', 'MET', 'LY6E', 'FOLH1', 'ERBB2', 'FAP', 'DLL3', 'CLDN6', 'SSTR2', 'EGFR']
kc = ['#E74C3C', '#27AE60', '#2980B9', '#1ABC9C', '#F39C12', '#8E44AD', '#E67E22', '#16A085', '#C0392B', '#D35400']
rx = ax[0].get_xlim()[1] * 0.82
for i in range(len(kn)):
    r = t[t['gene'].str.upper() == kn[i].upper()]
    if len(r) == 0:
        continue
    tv = r['tsi'].values[0]
    ax[0].axvline(x=tv, color=kc[i], linewidth=1.8, linestyle='--', alpha=0.8)
    yp = 0.94 - i * 0.09
    ax[0].text(rx, yp, kn[i], fontsize=7, color=kc[i], fontweight='bold', ha='left', va='center', transform=ax[0].get_xaxis_transform())
    ax[0].plot([tv, rx], [yp, yp], color=kc[i], alpha=0.2, linewidth=0.5, transform=ax[0].get_xaxis_transform())
ax[0].set_xlabel('TSI')
ax[0].set_ylabel('Count')
dm = e[e['depmap_score'].notna()]
sc = ax[1].scatter(dm['tsi_norm'], dm['depmap_score'], c=dm['ldt_norm'], cmap='RdYlGn', alpha=0.5, s=15, edgecolors='none')
ax[1].axhline(y=-0.5, color=RE, linestyle='--', linewidth=1, alpha=0.6)
ax[1].set_xlabel('TSI (norm)')
ax[1].set_ylabel('DepMap Chronos')
plt.tight_layout(pad=1.5)
save(fig3, 'figure3_validation.eps')
plt.close()

# Fig 4
print('Fig4...')
fig4, ax = plt.subplots(1, 2, figsize=(13, 5))
t30 = e.head(30)[::-1]
ec = 'tumor_median_log2' if 'tumor_median_log2' in t30.columns else 'tsi_norm'
ax[0].barh(range(30), t30[ec], height=0.65, color=BL, alpha=0.8)
ax[0].set_yticks(range(30))
ax[0].set_yticklabels(t30['gene'].values, fontsize=7)
ax[0].set_xlabel('Tumor log2(TPM)')
t15 = e.head(15)['gene'].tolist()
t15_ct = [g for g in t15 if g in ct['gene'].values]
cs = ct[ct['gene'].isin(t15_ct)]
ac = ['Colorectal', 'Lung Adeno', 'Pancreatic', 'Prostate', 'Thyroid',
      'Breast', 'Bladder', 'Kidney Clear Cell', 'Kidney Papillary',
      'Ovarian', 'Stomach', 'Liver', 'Head&Neck', 'Melanoma', 'Uterine',
      'Esophageal', 'Cervical', 'Sarcoma', 'Glioblastoma', 'Lung Squamous']
hd = pd.DataFrame(index=t15_ct, columns=ac)
for _, rw in cs.iterrows():
    g = rw['gene']
    if g in hd.index:
        if rw['best_cancer'] in hd.columns:
            hd.loc[g, rw['best_cancer']] = rw['best_log2tpm']
        if rw['second_cancer'] in hd.columns:
            hd.loc[g, rw['second_cancer']] = rw['second_log2tpm']
hd = hd.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
sns.heatmap(hd, annot=True, fmt='.1f', cmap='YlOrRd', linewidths=0.3,
    ax=ax[1], cbar_kws={'label': 'log2 TPM', 'shrink': 0.7}, annot_kws={'fontsize': 6})
ax[1].set_xlabel('Cancer Type', fontsize=8)
ax[1].set_ylabel('Gene', fontsize=8)
ax[1].tick_params(labelsize=7)
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=45, ha='right')
plt.tight_layout(pad=1.5)
save(fig4, 'figure4_expression_features.eps')
plt.close()

# Fig 5
print('Fig5...')
fig5, ax = plt.subplots(1, 2, figsize=(11, 4.5))
v = e.dropna(subset=['mean_plddt', 'ldt_norm'])
for rt in ['LYS', 'CYS', 'TYR', 'SER']:
    s = v[v['top_res_type'] == rt]
    ax[0].scatter(s['mean_plddt'], s['ldt_norm'], c=RC[rt], label=rt, alpha=0.4, s=10, edgecolors='none')
ax[0].set_xlabel('Mean pLDDT')
ax[0].set_ylabel('LDT Score (norm)')
ax[0].legend(frameon=False, fontsize=7, ncol=2)
ax[1].hist(f['n_pockets'], bins=60, color=PU, alpha=0.8, edgecolor='white', linewidth=0.3)
ax[1].set_xlabel('Pockets')
ax[1].set_ylabel('Count')
plt.tight_layout(pad=1.5)
save(fig5, 'figure5_quality.eps')
plt.close()

# Fig 6 Cancer heatmap (reuse hd)
print('Fig6...')
fig6, ax6 = plt.subplots(figsize=(14, 6))
sns.heatmap(hd, annot=True, fmt='.1f', cmap='YlOrRd', linewidths=0.3, ax=ax6,
    cbar_kws={'label': 'log2 TPM'}, annot_kws={'fontsize': 6})
ax6.set_title('Cancer-Type Specific Expression')
ax6.set_xlabel('Cancer Type')
ax6.set_ylabel('Gene')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
save(fig6, 'figure6_cancer_heatmap.eps')
plt.close()

# Fig 7 Benchmark
print('Fig7...')
fig7, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
ks = {'TACSTD2', 'MET', 'LY6E', 'FAP', 'FOLH1', 'ERBB2', 'EGFR', 'DLL3',
      'CLDN6', 'SSTR2', 'CD46', 'LRRC15', 'NECTIN4', 'STEAP1', 'CA9', 'CLDN18', 'MSLN'}
for name, col, color in [('TSI-only', 'tsi_norm', '#2471A3'),
                          ('LDT-only', 'ldt_norm', '#E67E22'),
                          ('Full model', 'enhanced_score', '#E74C3C')]:
    pos = {g for g in ks if g.upper() in {x.upper() for x in e['gene']}}
    labs = np.array([1 if r['gene'].upper() in pos else 0 for _, r in e.iterrows()])
    scs = e[col].values
    fpr, tpr, _ = roc_curve(labs, scs)
    auc = roc_auc_score(labs, scs)
    prec, rec, _ = precision_recall_curve(labs, scs)
    ap = average_precision_score(labs, scs)
    a1.plot(fpr, tpr, lw=2, color=color, label=name + ' AUC=' + str(round(auc, 2)))
    a2.plot(rec, prec, lw=2, color=color, label=name + ' AP=' + str(round(ap, 3)))
a1.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.3)
a1.set_xlabel('FPR')
a1.set_ylabel('TPR')
a1.set_title('ROC')
a1.legend(fontsize=7, frameon=False)
bl = 14 / len(e)
a2.axhline(y=bl, color='k', linestyle='--', lw=0.8, alpha=0.3)
a2.set_xlabel('Recall')
a2.set_ylabel('Precision')
a2.set_title('PR')
a2.legend(fontsize=7, frameon=False)
plt.tight_layout(pad=1.5)
save(fig7, 'figure_benchmark_roc.eps')
plt.close()

# Fig 8 Ablation
print('Fig8...')
fig8, ax8 = plt.subplots(figsize=(7, 3.5))
data = {'-TSI': -0.125, '-Pocket': -0.024, '-pLDDT': -0.002, '-LDT': 0.139}
comps = list(data.keys())[::-1]
vals = [data[c] for c in comps]
clrs = ['#3498DB' if v < 0 else '#E74C3C' for v in vals]
bars = ax8.barh(comps, vals, color=clrs, height=0.5, alpha=0.85)
ax8.axvline(x=0, color='black', linewidth=0.8)
ax8.set_xlabel('AUROC(without) - AUROC(full model)')
for bar, v, c in zip(bars, vals, clrs):
    if abs(v) > 0.1:
        ax8.text(v / 2, bar.get_y() + bar.get_height() / 2, str(round(v, 3)),
            va='center', ha='center', fontsize=8, fontweight='bold', color='white')
    else:
        ax8.text(v - 0.005 if v < 0 else v + 0.005,
            bar.get_y() + bar.get_height() / 2, str(round(v, 3)),
            va='center', ha='right' if v < 0 else 'left',
            fontsize=8, fontweight='bold', color=c)
ax8.spines['top'].set_visible(False)
ax8.spines['right'].set_visible(False)
from matplotlib.patches import Patch
ax8.legend(handles=[Patch(color='#E74C3C', label='Positive = removal improves'),
    Patch(color='#3498DB', label='Negative = removal degrades')], fontsize=7, frameon=False)
plt.tight_layout()
save(fig8, 'figure_ablation.eps')
plt.close()

# Fig 9 Filter cascade
print('Fig9...')
fig9, (a1, a2) = plt.subplots(2, 1, figsize=(10, 5), height_ratios=[2, 1])
sa = ['Surface\nProteome', 'Expression\nData', 'AlphaFold\nStructure',
      'Pockets\nDetected', 'LDT\nNucleophiles', 'Ligandable\n(>100 A3)']
ca = [2799, 2663, 2490, 2473, 2378, 2159]
clrs_a = ['#3498DB', '#2980B9', '#1ABC9C', '#27AE60', '#2ECC71', '#F39C12']
for i in range(len(sa)):
    a1.bar(i, ca[i], color=clrs_a[i], width=0.6, edgecolor='white')
    a1.text(i, ca[i] + 30, str(ca[i]), ha='center', fontsize=8, fontweight='bold')
a1.set_xticks(range(len(sa)))
a1.set_xticklabels(sa, fontsize=7)
a1.set_ylabel('Proteins')
a1.set_title('A. Global structural-chemical analysis')
sb = ['Top 200 Ranked', 'Extracellular Confirmed']
cb = [200, 157]
clrs_b = ['#E74C3C', '#27AE60']
for i in range(len(sb)):
    a2.bar(i, cb[i], color=clrs_b[i], width=0.4, edgecolor='white')
    a2.text(i, cb[i] + 5, str(cb[i]), ha='center', fontsize=9, fontweight='bold')
a2.set_xticks(range(len(sb)))
a2.set_xticklabels(sb, fontsize=8)
a2.set_ylabel('Proteins')
a2.set_title('B. Post hoc extracellular validation (UniProt topology)')
plt.tight_layout(pad=1.5)
save(fig9, 'figure_filter_cascade.eps')
plt.close()

# Fig 10 Tool comparison
print('Fig10...')
mu = pd.read_csv(BASE / 'results/tcga_mutation_frequencies.csv')
co = pd.read_csv(BASE / 'results/cross_species_conservation.csv')
ot = pd.read_csv(BASE / 'results/opentargets_disease.csv')
fig10, axx = plt.subplots(1, 2, figsize=(13, 4.5))
tools = ['LDT-TargetDB', 'DrugMap', 'ImmunoTar', 'TCSA']
features = [19, 11, 9, 8]
bar_colors = [RE, BL, GR, '#7F8C8D']
bars = axx[0].bar(tools, features, color=bar_colors, width=0.55, edgecolor='white', linewidth=0.5)
for b, c in zip(bars, features):
    axx[0].text(b.get_x() + b.get_width() / 2, b.get_height() + 0.3,
        str(c), ha='center', fontweight='bold', fontsize=13)
axx[0].set_ylabel('Features (of 20)')
axx[0].set_ylim(0, 23)
cats = ['TSI High\n(>0.7)', 'LDT Score\n(>0.5)', 'Low Mutation\n(<30)',
        'Conserved\n(>5 models)', 'Cancer-Assoc\n(OpenTargets)']
vals = [(e['tsi_norm'] > 0.7).sum(), (e['ldt_norm'] > 0.5).sum(),
    (mu['total_mutations'] < 30).sum() if 'total_mutations' in mu.columns else 0,
    (co['conservation_level'] == 'High (>5)').sum() if 'conservation_level' in co.columns else 0,
    (ot['n_cancer_associations'] > 0).sum() if 'n_cancer_associations' in ot.columns else 0]
axx[1].barh(cats, vals, color=RE, height=0.55, alpha=0.85)
for i in range(len(cats)):
    axx[1].text(vals[i] + 2, i, str(vals[i]), va='center', fontweight='bold', fontsize=10)
axx[1].set_xlabel('Qualifying Targets')
axx[1].set_xlim(0, max(vals) * 1.35)
plt.tight_layout(pad=1.5)
save(fig10, 'figure8_tool_comparison.eps')
plt.close()

# Fig 11 Workflow (user PNG -> EPS)
print('Fig11...')
wf = Image.open(str(BASE / 'results/figures/LDT_TargetDB_research_pipeline_900DPI.png')).convert('RGB')
wf.save(str(FIG / 'LDT_TargetDB_research_pipeline_900DPI.eps'), 'EPS', resolution=600)
print('  LDT_TargetDB_research_pipeline_900DPI.eps saved')

print('\n=== All 11 EPS files saved ===')
for f in sorted(FIG.glob('*.eps')):
    print('  {}: {:.1f} MB'.format(f.name, f.stat().st_size / 1024 / 1024))
