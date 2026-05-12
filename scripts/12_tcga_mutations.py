"""
TCGA突变频率分析：FireBrowse API → 33种癌种 × TOP20靶点
"""
import pandas as pd
import numpy as np
import requests
import time
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
RES = BASE / "results"

# TCGA cohorts (33 cancer types)
# FireBrowse cohort codes
TCGA_COHORTS = [
    "ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA",
    "GBM","GBMLGG","HNSC","KICH","KIPAN","KIRC","KIRP","LAML","LGG",
    "LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
    "SARC","SKCM","STAD","STES","TGCT","THCA","THYM","UCEC","UCS","UVM",
]

# 加载TOP靶点
enhanced = pd.read_csv(RES / "enhanced_final_ranking.csv")
top_genes = enhanced.head(20)["gene"].tolist()

print(f"Querying {len(top_genes)} genes across {len(TCGA_COHORTS)} TCGA cohorts...")

# ============================================================
# 查询每个基因×癌种的突变频率
# ============================================================
mutation_data = []

for i, gene in enumerate(top_genes):
    gene_mutations = {"gene": gene}

    for cohort in TCGA_COHORTS:
        try:
            url = "http://firebrowse.org/api/v1/Analyses/Mutation/MAF"
            params = {"cohort": cohort, "gene": gene, "page_size": 1}
            r = requests.get(url, params=params, timeout=10)

            if r.status_code == 200:
                data = r.json()
                maf = data.get("MAF", [])
                gene_mutations[cohort] = len(maf)
            else:
                gene_mutations[cohort] = -1
        except:
            gene_mutations[cohort] = -1

        time.sleep(0.05)  # 轻量限速

    # 汇总统计
    vals = [v for v in gene_mutations.values() if isinstance(v, (int, float)) and v >= 0]
    total_mutations = sum(vals)
    n_cohorts_with_mutations = sum(1 for v in vals if v > 0)
    max_mutations = max(vals) if vals else 0

    mutation_data.append({
        "gene": gene,
        "total_mutations": total_mutations,
        "n_cohorts_mutated": n_cohorts_with_mutations,
        "n_cohorts_queried": len(TCGA_COHORTS),
        "max_per_cohort": max_mutations,
        "mutation_frequency": round(total_mutations / len(TCGA_COHORTS), 1),
    })

    if (i + 1) % 5 == 0:
        print(f"  {i+1}/{len(top_genes)} done: {gene} mutations={total_mutations}")

mut_df = pd.DataFrame(mutation_data)
mut_df = mut_df.sort_values("total_mutations", ascending=False)

print(f"\n=== TCGA Mutation Frequency Results ===")
print(f"{'Gene':<15} {'Total Mutations':<18} {'Cohorts Mutated':<16} {'Max/Cohort':<12} {'Assessment'}")
print("-" * 75)
for _, row in mut_df.iterrows():
    if row["total_mutations"] < 20:
        assessment = "LOW — very stable target"
    elif row["total_mutations"] < 100:
        assessment = "MODERATE — acceptable"
    else:
        assessment = "HIGH — may affect targeting"

    print(f"{row['gene']:<15} {row['total_mutations']:<18} "
          f"{row['n_cohorts_mutated']}/{row['n_cohorts_queried']:<12} "
          f"{row['max_per_cohort']:<12} {assessment}")

# 标注口袋残基是否受影响
print(f"\n=== Pocket Residue Mutation Risk ===")
for _, row in enhanced.head(20).iterrows():
    gene = row["gene"]
    top_res = row.get("top_res_type", "")
    mr = mut_df[mut_df["gene"] == gene]
    if len(mr) > 0:
        total = mr["total_mutations"].values[0]
        risk = "Low" if total < 50 else "Moderate"
        print(f"  {gene:15s} Pocket: {top_res:4s} Mutations: {total:5d} Risk: {risk}")

mut_df.to_csv(RES / "tcga_mutation_frequencies.csv", index=False)
print(f"\nSaved to: {RES}/tcga_mutation_frequencies.csv")
print("Done.")
