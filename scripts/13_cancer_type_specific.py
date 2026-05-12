"""
癌种特异性分析：利用cBioPortal样本映射 + 本地TCGA表达数据
"""
import pandas as pd
import numpy as np
import gzip, json
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")

# 加载样本→癌种映射
with open(BASE / "data/expression/tcga_cancer_type_map.json") as f:
    cancer_samples = json.load(f)

# 加载基因映射
with open(BASE / "data/expression/ensg2symbol_cache.json") as f:
    ensg2symbol = json.load(f)
symbol2ensg = {v: k for k, v in ensg2symbol.items()}

# 加载增强排名TOP30
enhanced = pd.read_csv(BASE / "results/enhanced_final_ranking.csv")
top_genes = enhanced.head(30)["gene"].tolist()

# 构建样本→癌种的快速查找
sample_to_ct = {}
for ct, samples in cancer_samples.items():
    for s in samples:
        sample_to_ct[s] = ct

# 读取TCGA header获取样本顺序
with gzip.open(BASE / "data/expression/tcga_RSEM_gene_tpm.gz", "rt") as f:
    header = f.readline().strip().split("\t")[1:]

# 构建样本索引→癌种
sample_ct_list = [sample_to_ct.get(s, "Unknown") for s in header]
ct_counts = pd.Series(sample_ct_list).value_counts()
major_cts = [ct for ct in ct_counts[ct_counts >= 20].index if ct != "Unknown"]
print(f"Cancer types with >=20 samples: {len(major_cts)}")

# 逐行提取目标基因表达
target_ensgs = set()
for g in top_genes:
    if g in symbol2ensg:
        target_ensgs.add(symbol2ensg[g])

gene_expr = {}
with gzip.open(BASE / "data/expression/tcga_RSEM_gene_tpm.gz", "rt") as f:
    f.readline()
    for line in f:
        parts = line.strip().split("\t")
        ensg = parts[0].split(".")[0]
        if ensg in target_ensgs:
            vals = np.array([float(x) if x else np.nan for x in parts[1:]])
            gene_expr[ensg] = vals
            if len(gene_expr) >= len(target_ensgs):
                break

print(f"Genes extracted: {len(gene_expr)}")

# 计算每个基因×癌种的中位表达
results = []
for gene in top_genes:
    ensg = symbol2ensg.get(gene)
    if not ensg or ensg not in gene_expr:
        continue

    vals = gene_expr[ensg]
    ct_medians = {}
    for ct in major_cts:
        ct_idx = [i for i, c in enumerate(sample_ct_list) if c == ct]
        ct_vals = vals[ct_idx]
        ct_vals = ct_vals[~np.isnan(ct_vals)]
        if len(ct_vals) >= 10:
            ct_medians[ct] = np.median(ct_vals)

    if not ct_medians:
        continue

    sorted_ct = sorted(ct_medians.items(), key=lambda x: -x[1])
    top3 = sorted_ct[:3]
    n_positive = sum(1 for v in ct_medians.values() if v > 0)

    results.append({
        "gene": gene,
        "best_cancer": top3[0][0],
        "best_log2tpm": round(top3[0][1], 1),
        "second_cancer": top3[1][0] if len(top3) > 1 else "",
        "second_log2tpm": round(top3[1][1], 1) if len(top3) > 1 else 0,
        "third_cancer": top3[2][0] if len(top3) > 2 else "",
        "third_log2tpm": round(top3[2][1], 1) if len(top3) > 2 else 0,
        "n_cancers_positive": n_positive,
        "n_cancers_tested": len(ct_medians),
        "expression_range": round(sorted_ct[0][1] - sorted_ct[-1][1], 1),
    })

ct_df = pd.DataFrame(results).sort_values("best_log2tpm", ascending=False)

print(f"\n=== Cancer-Type Specific Target Expression ===\n")
print(f"{'Gene':<14} {'Best Cancer':<22} {'log2TPM':<8} {'2nd Cancer':<22} {'2nd':<8} {'#Cancers>0'}")
print("-" * 85)
for _, row in ct_df.iterrows():
    print(f"{row['gene']:<14} {row['best_cancer']:<22} {row['best_log2tpm']:<8.1f} "
          f"{row['second_cancer']:<22} {row['second_log2tpm']:<8.1f} "
          f"{row['n_cancers_positive']}/{row['n_cancers_tested']}")

# 统计每个癌种的最佳靶点
print(f"\n=== Best Targets Per Cancer Type ===\n")
for ct in sorted(major_cts):
    best_in_ct = ct_df[ct_df["best_cancer"] == ct]
    if len(best_in_ct) > 0:
        genes_str = ", ".join(best_in_ct["gene"].head(3).tolist())
        print(f"  {ct:25s}: {genes_str}")
    else:
        # 找这个癌种表达最高的靶点
        pass

ct_df.to_csv(BASE / "results/cancer_type_specific_expression.csv", index=False)
print(f"\nSaved to results/cancer_type_specific_expression.csv")
print("Done.")
