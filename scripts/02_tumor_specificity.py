"""
模块2：肿瘤特异性表达筛选 (内存优化版)
逐行读取TCGA/GTEx，只提取表面基因的表达值
"""
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import json

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
DATA = BASE / "data"

# ============================================================
# 1. ENSG → Gene Symbol 映射
# ============================================================
print("Building ENSG → Symbol mapping via mygene...")

mapping_cache = DATA / "expression" / "ensg2symbol_cache.json"
if mapping_cache.exists():
    print("  Loading cached mapping...")
    with open(mapping_cache) as f:
        ensg2symbol = json.load(f)
else:
    import mygene
    mg = mygene.MyGeneInfo()

    # 先收集所有ENSG ID
    all_ensgs = set()
    for fpath in [DATA / "expression" / "tcga_RSEM_gene_tpm.gz",
                  DATA / "expression" / "gtex_RSEM_gene_tpm.gz"]:
        with gzip.open(fpath, "rt") as f:
            f.readline()  # header
            for line in f:
                ensg = line.split("\t")[0].split(".")[0]
                all_ensgs.add(ensg)
                if len(all_ensgs) >= 70000:
                    break

    ensg_list = list(all_ensgs)
    print(f"  Querying {len(ensg_list)} ENSG IDs...")

    ensg2symbol = {}
    for i in range(0, len(ensg_list), 1000):
        batch = ensg_list[i:i+1000]
        try:
            results = mg.querymany(batch, scopes="ensembl.gene",
                                   fields="symbol", species="human")
            for r in results:
                sym = r.get("symbol")
                if sym and isinstance(sym, str):
                    ensg2symbol[r["query"]] = sym
        except Exception as e:
            print(f"    Error batch {i//1000}: {e}")
        if (i + 1000) % 10000 == 0:
            print(f"    {i+1000}/{len(ensg_list)}, mapped: {len(ensg2symbol)}")

    with open(mapping_cache, "w") as f:
        json.dump(ensg2symbol, f)

print(f"  Mapped: {len(ensg2symbol)} genes")
symbol2ensg = {v: k for k, v in ensg2symbol.items()}

# ============================================================
# 2. 加载表面基因 → 找对应的ENSG
# ============================================================
surface_genes_raw = pd.read_csv(DATA / "surfaceome" / "core_surface_genes.txt", header=None)[0].tolist()
valid_genes = []
for g in surface_genes_raw:
    g = str(g).strip()
    if g in symbol2ensg:
        valid_genes.append(g)

target_ensgs = set(symbol2ensg[g] for g in valid_genes)
print(f"\nSurface genes mappable: {len(valid_genes)}/{len(surface_genes_raw)}")

# ============================================================
# 3. 逐行读取，只提取目标基因
# ============================================================
def extract_target_genes(filepath, target_ensgs):
    """逐行读取gzip压缩的表达矩阵，返回 {ensg: [values], samples: [sample_ids]}"""
    with gzip.open(filepath, "rt") as f:
        header = f.readline().strip().split("\t")
        samples = header[1:]

        gene_data = {}
        for line in f:
            parts = line.strip().split("\t")
            ensg = parts[0].split(".")[0]
            if ensg in target_ensgs:
                vals = np.array([float(x) if x else np.nan for x in parts[1:]])
                gene_data[ensg] = vals
                if len(gene_data) >= len(target_ensgs):
                    break
    return samples, gene_data

print("Extracting TCGA data...")
tcga_samples, tcga_data = extract_target_genes(
    DATA / "expression" / "tcga_RSEM_gene_tpm.gz", target_ensgs)
print(f"  Samples: {len(tcga_samples)}, Genes extracted: {len(tcga_data)}")

print("Extracting GTEx data...")
gtex_samples, gtex_data = extract_target_genes(
    DATA / "expression" / "gtex_RSEM_gene_tpm.gz", target_ensgs)
print(f"  Samples: {len(gtex_samples)}, Genes extracted: {len(gtex_data)}")

# ============================================================
# 4. 分类样本
# ============================================================
def classify_tcga(sample_id):
    parts = sample_id.split("-")
    if len(parts) < 4:
        return "unknown", "unknown"
    cancer_type = parts[1]
    sample_code = parts[3]
    if sample_code.startswith("01"):
        return cancer_type, "tumor"
    elif sample_code.startswith("11"):
        return cancer_type, "normal"
    else:
        return cancer_type, "other"

tcga_tumor_idx = []
tcga_normal_idx = []
tcga_ct_map = {}

for i, s in enumerate(tcga_samples):
    ct, stype = classify_tcga(s)
    if stype == "tumor":
        tcga_tumor_idx.append(i)
        tcga_ct_map[i] = ct
    elif stype == "normal":
        tcga_normal_idx.append(i)

print(f"\nTCGA tumor: {len(tcga_tumor_idx)}, normal: {len(tcga_normal_idx)}")
print(f"GTEx samples: {len(gtex_samples)} (all normal tissue)")

# ============================================================
# 5. 计算TSI
# ============================================================
print("\n=== Calculating TSI ===")

results = []
for gene in valid_genes:
    ensg = symbol2ensg[gene]

    tcga_arr = tcga_data.get(ensg)
    gtex_arr = gtex_data.get(ensg)

    if tcga_arr is None:
        continue

    # TCGA肿瘤值
    tumor_vals = tcga_arr[tcga_tumor_idx] if tcga_tumor_idx else np.array([])
    tumor_vals = tumor_vals[~np.isnan(tumor_vals)]

    # TCGA正常值
    normal_vals = tcga_arr[tcga_normal_idx] if tcga_normal_idx else np.array([])
    normal_vals = normal_vals[~np.isnan(normal_vals)]

    # GTEx正常值
    gtex_vals = np.array([])
    if gtex_arr is not None:
        gtex_vals = gtex_arr[~np.isnan(gtex_arr)]

    all_normal = np.concatenate([normal_vals, gtex_vals])

    if len(tumor_vals) < 50:
        continue

    tumor_median = np.median(tumor_vals)
    normal_median = np.median(all_normal) if len(all_normal) > 0 else -10

    # FC (log2 scale)
    fc_log2 = tumor_median - normal_median

    # GTEx中高于背景的组织数 (-9.9658 = log2(0.001), 用-3 ≈ TPM~0.125)
    gtex_positive = np.sum(gtex_vals > -3) if len(gtex_vals) > 0 else 0

    # 肿瘤阳性率 (>0 ≈ TPM>1)
    tumor_positive_rate = np.mean(tumor_vals > 0)

    # TSI (pan-cancer)
    tsi = fc_log2 * tumor_positive_rate / np.sqrt(gtex_positive + 1)

    results.append({
        "gene": gene,
        "ensg": ensg,
        "tumor_median_log2": round(tumor_median, 3),
        "normal_median_log2": round(normal_median, 3),
        "fc_log2": round(fc_log2, 2),
        "gtex_positive_count": int(gtex_positive),
        "tumor_positive_rate": round(tumor_positive_rate, 3),
        "tsi": round(tsi, 3),
    })

# ============================================================
# 6. 排序 + 输出
# ============================================================
tsi_df = pd.DataFrame(results).sort_values("tsi", ascending=False).reset_index(drop=True)
tsi_df["tsi_norm"] = ((tsi_df["tsi"] - tsi_df["tsi"].min()) /
                       (tsi_df["tsi"].max() - tsi_df["tsi"].min() + 1e-10))
tsi_df["rank"] = range(1, len(tsi_df) + 1)

print(f"\nTSI calculated: {len(tsi_df)} genes")
print(f"TSI range: {tsi_df['tsi'].min():.2f} - {tsi_df['tsi'].max():.2f}\n")
print("Top 30:")
for _, row in tsi_df.head(30).iterrows():
    print(f"  {row['rank']:3d}. {row['gene']:15s} TSI={row['tsi']:7.2f}  "
          f"FC={row['fc_log2']:5.1f}  PosRate={row['tumor_positive_rate']:.2f}  "
          f"NormalCt={row['gtex_positive_count']:3d}")

# 保存
tsi_df.to_csv(DATA / "expression" / "tsi_results.csv", index=False)
tsi_df.to_pickle(DATA / "expression" / "tsi_results.pkl")

top200 = tsi_df.head(200)
top200["gene"].to_csv(DATA / "surfaceome" / "candidate_genes_top200.txt",
                      index=False, header=False)

print(f"\nDone. Top 200 genes → candidate_genes_top200.txt")
