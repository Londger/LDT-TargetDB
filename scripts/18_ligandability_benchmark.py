"""
Pocket ligandability upgrade + rigorous benchmarking
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
RES = BASE / "results"
FIG = BASE / "results/figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial'],'font.size':9,
                      'axes.titlesize':11,'axes.labelsize':10,'axes.spines.top':False,
                      'axes.spines.right':False,'savefig.dpi':600,'savefig.bbox':'tight',
                      'figure.facecolor':'white','axes.facecolor':'white'})

# Load data
enhanced = pd.read_csv(RES / "enhanced_final_ranking.csv")
full = pd.read_csv(RES / "full_proteome_ranking.csv")

# ============================================================
# 1. Ligandability Upgrade
# ============================================================
print("=== 1. Pocket Ligandability Upgrade ===")

# Apply geometric filters to full proteome data
# Filter: volume >= 100 Å³ (typical drug-like pocket minimum)
# Keep targets with at least one "ligandable" pocket
VOLUME_THRESHOLD = 100  # Å³

ligandable = full.copy()
ligandable["is_ligandable"] = ligandable["max_pocket_volume"] >= VOLUME_THRESHOLD
ligandable["pocket_volume_log10"] = np.log10(ligandable["max_pocket_volume"].clip(lower=1))

n_ligandable = ligandable["is_ligandable"].sum()
n_total = len(ligandable)
print(f"  Volume >= {VOLUME_THRESHOLD} Å³: {n_ligandable}/{n_total} ({n_ligandable/n_total*100:.1f}%)")
print(f"  Mean volume: {ligandable['max_pocket_volume'].mean():.0f} Å³, Median: {ligandable['max_pocket_volume'].median():.0f} Å³")

# New scoring: LDT × ligandability penalty
# Ligandability penalty: volume-based sigmoid
ligandable["ligandability_score"] = 1 / (1 + np.exp(-(ligandable["max_pocket_volume"] - VOLUME_THRESHOLD) / 50))

# Update LDT score with ligandability
ligandable["ldt_ligandable"] = ligandable["ldt_score"] * ligandable["ligandability_score"]

# Update composite score
ligandable["pocket_norm_v2"] = np.clip(ligandable["max_pocket_volume"] / ligandable["max_pocket_volume"].max(), 0, 1)
ligandable["ldt_norm_v2"] = np.clip(ligandable["ldt_ligandable"] / ligandable["ldt_ligandable"].max(), 0, 1)

# Merge TSI
tsi = pd.read_csv(BASE / "data/expression/tsi_results.csv")
tsi_lookup = dict(zip(tsi["gene"], tsi["tsi_norm"]))
ligandable["tsi_norm"] = ligandable["gene"].map(tsi_lookup).fillna(0)

ligandable["final_score_v3"] = (0.4 * ligandable["tsi_norm"].fillna(0) +
                                 0.3 * ligandable["pocket_norm_v2"].fillna(0) +
                                 0.3 * ligandable["ldt_norm_v2"].fillna(0))

ligandable = ligandable.sort_values("final_score_v3", ascending=False).reset_index(drop=True)
ligandable["rank_v3"] = range(1, len(ligandable)+1)

print(f"\n  Top 10 after ligandability upgrade:")
for _, row in ligandable.head(10).iterrows():
    print(f"  {row['rank_v3']:3d}. {row['gene']:15s} Score={row['final_score_v3']:.3f} "
          f"Vol={row['max_pocket_volume']:.0f}Å³ Lig={row['ligandability_score']:.2f}")

# ============================================================
# 2. Rigorous Benchmark
# ============================================================
print("\n=== 2. Benchmark with AUROC ===")

# Known positive set (17 nuclear medicine targets)
known_positives = {"TACSTD2","MET","LY6E","FAP","FOLH1","ERBB2","EGFR","DLL3","CLDN6",
                   "SSTR2","MSLN","CD46","LRRC15","NECTIN4","STEAP1","CA9","CLDN18"}

# Build labeled dataset
def build_labeled_dataset(ranked_df, score_col, positive_set):
    """Create binary labels and scores for AUC calculation"""
    genes = set(ranked_df["gene"].tolist())
    pos_set = {g for g in positive_set if g.upper() in {x.upper() for x in genes}}

    labels = []
    scores = []
    for _, row in ranked_df.iterrows():
        gene = row["gene"]
        labels.append(1 if gene.upper() in pos_set else 0)
        scores.append(row[score_col])
    return np.array(labels), np.array(scores)

from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

# Models to compare
models = {
    "Expression-only (TSI)": ("tsi_norm", enhanced),
    "LDT Chemistry-only": ("ldt_norm", enhanced),
    "Pocket volume-only": ("pocket_norm", enhanced),
    "Full model (original)": ("enhanced_score", enhanced),
    "Full model (v3 ligandable)": ("final_score_v3", ligandable),
}

benchmark_results = []
fig_roc, ax_roc = plt.subplots(figsize=(5,5))
fig_pr, ax_pr = plt.subplots(figsize=(5,5))

for model_name, (score_col, df) in models.items():
    labels, scores = build_labeled_dataset(df, score_col, known_positives)

    # Skip if no positives
    n_pos = labels.sum()
    if n_pos == 0:
        print(f"  {model_name}: 0 positives, skipping")
        continue

    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    baseline = n_pos / len(labels)

    # Top 1% enrichment
    top1_k = max(1, int(len(labels) * 0.01))
    top1_idx = np.argsort(scores)[-top1_k:]
    top1_hits = labels[top1_idx].sum()
    top1_enrich = top1_hits / (n_pos * 0.01) if n_pos > 0 else 0

    benchmark_results.append({
        "model": model_name,
        "AUROC": round(auroc, 3),
        "AUPRC": round(auprc, 3),
        "baseline": round(baseline, 3),
        "top1_enrich": round(top1_enrich, 1),
        "n_positives": int(n_pos),
        "n_total": len(labels),
    })

    # ROC curve
    fpr, tpr, _ = roc_curve(labels, scores)
    ax_roc.plot(fpr, tpr, lw=1.5, label=f"{model_name} (AUC={auroc:.2f})")

    # PR curve
    prec, rec, _ = precision_recall_curve(labels, scores)
    ax_pr.plot(rec, prec, lw=1.5, label=f"{model_name} (AP={auprc:.2f})")

    print(f"  {model_name}: AUROC={auroc:.3f}, AUPRC={auprc:.3f} (baseline={baseline:.3f}), "
          f"top1_enrich={top1_enrich:.1f}x, n_pos={int(n_pos)}")

# ROC styling
ax_roc.plot([0,1],[0,1],'k--',lw=0.8,alpha=0.3)
ax_roc.set_xlabel('False Positive Rate'), ax_roc.set_ylabel('True Positive Rate')
ax_roc.set_title('ROC: Known Nuclear Medicine Target Recovery')
ax_roc.legend(fontsize=7,frameon=False)

# PR styling
baseline_val = benchmark_results[0]["baseline"] if benchmark_results else 0.05
ax_pr.axhline(y=baseline_val, color='k', linestyle='--', lw=0.8, alpha=0.3, label=f'Baseline ({baseline_val:.3f})')
ax_pr.set_xlabel('Recall'), ax_pr.set_ylabel('Precision')
ax_pr.set_title('PR Curve: Known Nuclear Medicine Target Recovery')
ax_pr.legend(fontsize=7,frameon=False)

fig_roc.savefig(FIG / "figure_benchmark_roc.png")
fig_pr.savefig(FIG / "figure_benchmark_pr.png")
plt.close('all')

# ============================================================
# 3. Filter Cascade Plot
# ============================================================
print("\n=== 3. Filter Cascade ===")
fig_cascade, ax = plt.subplots(figsize=(8,4))

stages = ["Surface\nProteome","Expression\nData","AlphaFold\nStructure","Binding\nPockets","LDT\nNucleophiles","Ligandable\n(>100Å³)","Extracellular\nConfirmed"]
counts = [2799, 2663, 2490, 2473, 2378, n_ligandable, 157]
colors_cascade = ['#3498DB','#2980B9','#1ABC9C','#27AE60','#2ECC71','#F39C12','#E74C3C']

for i, (stage, count, color) in enumerate(zip(stages, counts, colors_cascade)):
    ax.bar(i, count, color=color, width=0.6, edgecolor='white', linewidth=0.5)
    ax.text(i, count + 30, str(count), ha='center', fontsize=8, fontweight='bold')
ax.set_xticks(range(len(stages)))
ax.set_xticklabels(stages, fontsize=7)
ax.set_ylabel('Number of Proteins')
ax.set_title('Filter Cascade: From Surface Proteome to Prioritized Targets')
ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)

fig_cascade.savefig(FIG / "figure_filter_cascade.png")
plt.close()

# ============================================================
# 4. Ablation Analysis
# ============================================================
print("\n=== 4. Ablation ===")

ablations = {
    "Full model": enhanced["enhanced_score"],
    "-TSI": 0.25*enhanced["pocket_norm"].fillna(0) + 0.25*enhanced["ldt_norm"].fillna(0) + 0.10*(1-enhanced["plddt_below_70"].fillna(100)/100) + 0.05*np.where(enhanced["depmap_label"]=="Essential",1,np.where(enhanced["depmap_label"]=="No data",0,0.3)),
    "-LDT": enhanced["tsi_norm"]*0.5 + enhanced["pocket_norm"]*0.2 + 0.10*(1-enhanced["plddt_below_70"].fillna(100)/100),
    "-Pocket": enhanced["tsi_norm"]*0.5 + enhanced["ldt_norm"]*0.3 + 0.10*(1-enhanced["plddt_below_70"].fillna(100)/100),
    "-pLDDT": enhanced["tsi_norm"]*0.4 + enhanced["pocket_norm"]*0.3 + enhanced["ldt_norm"]*0.3,
    "-DepMap": enhanced["tsi_norm"]*0.4 + enhanced["pocket_norm"]*0.3 + enhanced["ldt_norm"]*0.3,
}

ablation_results = []
for ab_name, scores in ablations.items():
    labels, scores_arr = build_labeled_dataset(
        pd.DataFrame({"gene": enhanced["gene"], "score": scores}), "score", known_positives)
    if labels.sum() > 0:
        auc = roc_auc_score(labels, scores_arr)
        ablation_results.append({"Ablation": ab_name, "AUROC": round(auc,3)})
        print(f"  {ab_name}: AUROC={auc:.3f}")

ablation_df = pd.DataFrame(ablation_results)

fig_ab, ax_ab = plt.subplots(figsize=(7,4))
full_auc = ablation_df[ablation_df["Ablation"]=="Full model"]["AUROC"].values[0]
auc_drops = []
for _, row in ablation_df.iterrows():
    if row["Ablation"] != "Full model":
        auc_drops.append({"component": row["Ablation"], "drop": full_auc - row["AUROC"]})

drop_df = pd.DataFrame(auc_drops).sort_values("drop", ascending=True)
ax_ab.barh(drop_df["component"], drop_df["drop"], color='#E74C3C', height=0.5, alpha=0.8)
ax_ab.set_xlabel('AUROC Drop When Removed')
ax_ab.set_title('Ablation Analysis: Impact of Removing Each Component')
for i, (_, row) in enumerate(drop_df.iterrows()):
    ax_ab.text(row["drop"]+0.005, i, f'{row["drop"]:.3f}', va='center', fontsize=8, fontweight='bold')

fig_ab.savefig(FIG / "figure_ablation.png")
plt.close()

# ============================================================
# 5. Save Results
# ============================================================
benchmark_df = pd.DataFrame(benchmark_results)
benchmark_df.to_csv(RES / "benchmark_results.csv", index=False)
ablation_df.to_csv(RES / "ablation_results.csv", index=False)
ligandable[["gene","rank_v3","final_score_v3","max_pocket_volume","ligandability_score"]].head(100).to_csv(
    RES / "ligandable_top100.csv", index=False)

print(f"\n=== Complete ===")
print(f"Generated: figure_benchmark_roc.png, figure_benchmark_pr.png, figure_filter_cascade.png, figure_ablation.png")
print(f"Saved: benchmark_results.csv, ablation_results.csv, ligandable_top100.csv")
