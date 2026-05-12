"""
模块5：增强分析 — 癌种特异性 + 已知靶点验证 + DepMap + pLDDT质控
"""
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import json

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
DATA = BASE / "data"
RES = BASE / "results"
RES.mkdir(exist_ok=True)

# ============================================================
# 1. 加载现有结果
# ============================================================
print("Loading results...")
tsi = pd.read_csv(DATA / "expression" / "tsi_results.csv")
final = pd.read_csv(BASE / "results" / "full_proteome_ranking.csv")
struct_map = pd.read_csv(DATA / "structures" / "all_structures_combined.csv")
failed = None  # 全量管道不追踪单个失败

# ============================================================
# 2. 已知核药靶点完整分析
# ============================================================
print("\n=== 2. Known Nuclear Medicine Targets ===")

known_targets = {
    "FOLH1": "PSMA", "FAP": "FAP", "SSTR2": "SSTR2", "SSTR5": "SSTR5",
    "GRPR": "GRPR", "DLL3": "DLL3", "STEAP1": "STEAP1", "CD46": "CD46",
    "LRRC15": "LRRC15", "LY6E": "LY6E", "NECTIN4": "Nectin-4",
    "TACSTD2": "Trop-2", "ERBB2": "HER2", "EGFR": "EGFR",
    "MET": "c-Met", "MSLN": "Mesothelin", "GPC3": "GPC3",
    "CA9": "CAIX", "CLDN6": "CLDN6", "CLDN18": "CLDN18.2",
}

# 统一大小写
tsi["gene_clean"] = tsi["gene"].str.strip()
final["gene_clean"] = final["gene"].str.strip()
known_list_upper = [k.upper() for k in known_targets]

# TSI的rank改名避免冲突
tsi.rename(columns={"rank": "tsi_rank"}, inplace=True)

known_df = tsi[tsi["gene"].str.upper().isin(known_list_upper)].copy()
final_cols = ["gene_clean", "rank", "final_score", "n_pockets", "n_nucleophiles",
              "ldt_score", "top_res_type"]
final_cols = [c for c in final_cols if c in final.columns]
known_df = known_df.merge(
    final[final_cols],
    left_on="gene_clean", right_on="gene_clean", how="left"
)
known_df["standard_name"] = known_df["gene"].str.upper().map(
    {k.upper(): v for k, v in known_targets.items()})
known_df = known_df.sort_values("tsi", ascending=False)

print(f"Known targets in data: {len(known_df)}/{len(known_targets)}")
in_final = known_df["rank"].notna().sum()
in_top50 = (known_df["rank"].fillna(999) <= 50).sum()
print(f"  With structural data: {in_final}/{len(known_df)} (in top200 TSI → had AlphaFold)")
print(f"  In LDT TOP 50: {in_top50}")
print(f"\n{'Gene':<12} {'Name':<14} {'TSI':<8} {'TSI Rank':<10} {'LDT Rank':<10} {'Score':<8} {'LDT':<10} {'Best Res'}")
print("-"*85)
for _, row in known_df.iterrows():
    r = f"{int(row['rank'])}" if pd.notna(row['rank']) else "N/A"
    tr = f"{int(row['tsi_rank'])}" if pd.notna(row['tsi_rank']) else "-"
    s = f"{row['final_score']:.3f}" if pd.notna(row['final_score']) else "-"
    l = f"{row['ldt_score']:.4f}" if pd.notna(row['ldt_score']) else "-"
    res = f"{row['top_res_type']}" if pd.notna(row.get('top_res_type')) else "-"
    print(f"{row['gene']:<12} {row['standard_name']:<14} {row['tsi']:.2f}    {tr:<10} {r:<10} {s:<8} {l:<10} {res}")

# ============================================================
# 3. 肿瘤表达特征（泛癌）
# ============================================================
print("\n=== 3. Pan-Cancer Expression Features ===")
# 使用TSI中的tumor_median_log2作为泛癌表达水平
tsi_expr = tsi[["gene", "tumor_median_log2", "tumor_positive_rate"]].copy()
print(f"Genes with tumor expression data: {len(tsi_expr)}")
print(f"Mean tumor log2(TPM): {tsi_expr['tumor_median_log2'].mean():.1f}")
print(f"Mean positive rate: {tsi_expr['tumor_positive_rate'].mean():.2f}")

# 为增强排名准备简化版表达特征
best_ct_df = pd.DataFrame({
    "gene": tsi_expr["gene"],
    "tumor_median_log2": tsi_expr["tumor_median_log2"],
    "tumor_positive_rate": tsi_expr["tumor_positive_rate"],
})

# ============================================================
# 4. DepMap CRISPR 必需性
# ============================================================
print("\n=== 4. DepMap CRISPR Essentiality ===")

depmap_raw = pd.read_csv(DATA / "expression" / "CRISPRGeneEffect.csv", index_col=0)

# DepMap is transposed: rows=cell lines, cols=genes. Transpose.
depmap = depmap_raw.T
# Column names: "A1BG (1)" -> extract "A1BG"
depmap.index = depmap.index.astype(str).str.extract(r"^(.+?)\s*\(")[0].fillna(
    pd.Series(depmap.index.astype(str), index=depmap.index))

print(f"DepMap: {depmap.shape[0]} genes x {depmap.shape[1]} cell lines")

depmap_index_clean = depmap.index.astype(str).str.strip()
depmap_gene_map = dict(zip(depmap_index_clean, depmap.index))

# 计算每个基因的中位必需性分数
gene_median_ess = {}
for clean_name, orig_name in depmap_gene_map.items():
    row = depmap.loc[orig_name]
    gene_median_ess[clean_name] = row.median()

# 跟我们的靶点交叉
top_gene_ess = {}
for gene in final["gene"]:
    ess = gene_median_ess.get(gene)
    if ess is not None:
        top_gene_ess[gene] = ess

# 标记必需基因 (Chronos score < -0.5 被认为必需)
essential_genes = [g for g, v in top_gene_ess.items() if v < -0.5]
print(f"Genes with DepMap data: {len(top_gene_ess)}")
print(f"Essential genes (score < -0.5): {len(essential_genes)}")

# TOP30 的DepMap状态
ess_status = []
for _, row in final.head(30).iterrows():
    ess = top_gene_ess.get(row["gene"])
    if ess is not None:
        label = "Essential" if ess < -0.5 else ("Context-dependent" if ess < -0.3 else "Non-essential")
    else:
        ess = np.nan
        label = "No data"
    ess_status.append({"gene": row["gene"], "depmap_median": ess, "depmap_label": label})

ess_df = pd.DataFrame(ess_status)
print(f"\nTOP30 DepMap status:")
for _, row in ess_df.iterrows():
    print(f"  {row['gene']:15s} DepMap={row['depmap_median']:.3f} ({row['depmap_label']})")

# ============================================================
# 5. pLDDT置信度评估
# ============================================================
print("\n=== 5. pLDDT Confidence ===")

# 解析AlphaFold PDB中的pLDDT (存于B-factor列)
plddt_scores = {}
for _, row in struct_map.iterrows():
    if pd.isna(row["pdb_file"]):
        continue
    pdb_file = Path(row["pdb_file"])
    if not pdb_file.exists():
        continue
    gene = row["gene"]
    bfactors = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                b = float(line[60:66].strip())
                bfactors.append(b)
    if bfactors:
        plddt_scores[gene] = {
            "mean_plddt": round(np.mean(bfactors), 1),
            "plddt_below_70": round(100 * sum(1 for b in bfactors if b < 70) / len(bfactors), 1),
            "plddt_below_50": round(100 * sum(1 for b in bfactors if b < 50) / len(bfactors), 1),
        }

plddt_df = pd.DataFrame(plddt_scores).T
plddt_df.index.name = "gene"
plddt_df = plddt_df.reset_index()

print(f"Structures with pLDDT data: {len(plddt_df)}")
print(f"Mean pLDDT distribution:")
print(f"  High (>90): {(plddt_df['mean_plddt'] > 90).sum()}")
print(f"  Medium (70-90): {((plddt_df['mean_plddt'] >= 70) & (plddt_df['mean_plddt'] <= 90)).sum()}")
print(f"  Low (<70): {(plddt_df['mean_plddt'] < 70).sum()}")

# ============================================================
# 6. 综合增强排名
# ============================================================
print("\n=== 6. Enhanced Final Ranking ===")

# 合并所有信息
enhanced_cols = ["gene", "rank", "final_score", "tsi_norm", "pocket_norm",
                 "ldt_norm", "ldt_score", "n_pockets", "n_nucleophiles",
                 "top_res_type", "uniprot_id"]
enhanced_cols = [c for c in enhanced_cols if c in final.columns]
enhanced = final[enhanced_cols].copy()

# +泛癌表达
enhanced = enhanced.merge(best_ct_df, on="gene", how="left")

# +DepMap
enhanced["depmap_score"] = enhanced["gene"].map(top_gene_ess)
enhanced["depmap_label"] = enhanced["depmap_score"].apply(
    lambda x: "Essential" if pd.notna(x) and x < -0.5 else
              ("Context-dependent" if pd.notna(x) and x < -0.3 else
               ("Non-essential" if pd.notna(x) else "No data")))

# +pLDDT
enhanced = enhanced.merge(plddt_df, on="gene", how="left")

enhanced["mean_plddt"] = enhanced["mean_plddt"].fillna(0)
enhanced["plddt_below_70"] = enhanced["plddt_below_70"].fillna(100)

# 增强综合得分
enhanced["enhanced_score"] = (
    0.35 * enhanced["tsi_norm"] +
    0.25 * enhanced["pocket_norm"] +
    0.25 * enhanced["ldt_norm"] +
    0.10 * (1 - enhanced["plddt_below_70"] / 100) +  # pLDDT置信度
    0.05 * np.where(enhanced["depmap_label"] == "Essential", 1,
             np.where(enhanced["depmap_label"] == "Context-dependent", 0.7,
             np.where(enhanced["depmap_label"] == "No data", 0, 0.3)))  # 必需性
)

enhanced = enhanced.sort_values("enhanced_score", ascending=False).reset_index(drop=True)
enhanced["enhanced_rank"] = range(1, len(enhanced) + 1)

print(f"\n{'Rank':<5} {'Gene':<14} {'Enhanced':<8} {'Original':<8} {'TSI':<7} {'LDT':<7} "
      f"{'DepMap':<16} {'pLDDT':<8} {'BestNuc'}")
print("-"*95)
for _, row in enhanced.head(30).iterrows():
    dm = f"{row['depmap_label']} ({row['depmap_score']:.2f})" if pd.notna(row['depmap_score']) else "No data"
    plddt_str = f"{row['mean_plddt']:.0f}" if pd.notna(row['mean_plddt']) else "-"
    nuc = f"{row['top_res_type']}" if pd.notna(row['top_res_type']) else "-"
    print(f"{row['enhanced_rank']:<5} {row['gene']:<14} {row['enhanced_score']:.3f}    "
          f"{row['final_score']:.3f}    {row['tsi_norm']:.3f}   {row['ldt_norm']:.3f}   "
          f"{dm:<16} {plddt_str:<8} {nuc}")

# ============================================================
# 7. 保存全部增强结果
# ============================================================
enhanced.to_csv(RES / "enhanced_final_ranking.csv", index=False)
enhanced.head(20).to_csv(RES / "enhanced_top20_targets.csv", index=False)
known_df.to_csv(RES / "known_targets_validation.csv", index=False)
ess_df.to_csv(RES / "depmap_essentiality.csv", index=False)
plddt_df.to_csv(RES / "plddt_confidence.csv", index=False)

print(f"\nAll results saved to: {RES}")
print("Done.")
