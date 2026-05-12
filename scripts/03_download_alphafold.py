"""
模块3前置：批量从AlphaFold DB下载候选靶点蛋白结构
使用UniProt ID映射 → AlphaFold PDB文件
"""
import pandas as pd
import time
import requests
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
STRUCT_DIR = BASE / "data" / "structures"
STRUCT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. 读取候选基因列表 + 获取UniProt ID映射
# ============================================================
# 加载TSI结果
tsi_df = pd.read_csv(BASE / "data" / "expression" / "tsi_results.csv")
print(f"Total genes with TSI: {len(tsi_df)}")

# 取TSI前200的基因做结构分析（控制下载量和计算量）
top200 = tsi_df.head(200)
print(f"Top 200 genes for structure analysis")

# 加载SURFY数据获取UniProt映射
surfy = pd.read_excel(
    BASE / "data" / "surfaceome" / "SURFY_surfaceome.xlsx",
    sheet_name="in silico surfaceome only",
    skiprows=1
)
surfy = surfy.rename(columns={"UniProt gene": "gene", "UniProt accession": "uniprot_id"})
surfy["uniprot_id"] = surfy["uniprot_id"].str.extract(r"([A-Z][0-9][A-Z0-9]{3,})")[0]
surfy = surfy.dropna(subset=["gene", "uniprot_id"])
gene_to_uniprot = dict(zip(surfy["gene"], surfy["uniprot_id"]))

# 获取候选基因的UniProt ID
candidates = []
for gene in top200["gene"]:
    uniprot = gene_to_uniprot.get(gene)
    if uniprot:
        candidates.append({"gene": gene, "uniprot_id": uniprot})
    else:
        print(f"  No UniProt ID for {gene}")

candidates_df = pd.DataFrame(candidates)
print(f"Genes with UniProt mapping: {len(candidates_df)}")

# ============================================================
# 2. 批量下载AlphaFold结构
# ============================================================
BASE_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.pdb"

downloaded = []
failed = []

for i, row in candidates_df.iterrows():
    uniprot = row["uniprot_id"]
    gene = row["gene"]
    output_file = STRUCT_DIR / f"AF-{uniprot}-F1-model_v6.pdb"

    if output_file.exists() and output_file.stat().st_size > 1000:
        downloaded.append({"gene": gene, "uniprot_id": uniprot, "pdb_file": str(output_file)})
        continue

    url = BASE_URL.format(uniprot=uniprot)
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and (r.text.startswith("HEADER") or r.text.startswith("ATOM")):
            output_file.write_text(r.text)
            downloaded.append({"gene": gene, "uniprot_id": uniprot, "pdb_file": str(output_file)})
        else:
            failed.append({"gene": gene, "uniprot_id": uniprot, "error": f"HTTP {r.status_code}"})
    except Exception as e:
        failed.append({"gene": gene, "uniprot_id": uniprot, "error": str(e)})

    # 进度
    if (i + 1) % 20 == 0:
        print(f"  Progress: {i+1}/{len(candidates_df)} (downloaded: {len(downloaded)}, failed: {len(failed)})")

    # AlphaFold DB 限速 (10 req/s)
    time.sleep(0.1)
    import sys; sys.stdout.flush()

# ============================================================
# 3. 保存结果
# ============================================================
print(f"\n=== Download Summary ===")
print(f"Downloaded: {len(downloaded)}")
print(f"Failed: {len(failed)}")

pd.DataFrame(downloaded).to_csv(STRUCT_DIR / "downloaded_structures.csv", index=False)
if failed:
    pd.DataFrame(failed).to_csv(STRUCT_DIR / "failed_downloads.csv", index=False)
    print("Failed genes (check failed_downloads.csv):")
    for f in failed[:10]:
        print(f"  {f['gene']}: {f.get('error', 'not found')}")

# 输出用于后续分析的基因-结构映射
with open(STRUCT_DIR / "structure_gene_map.txt", "w") as f:
    f.write("gene\tuniprot_id\tpdb_file\n")
    for d in downloaded:
        f.write(f"{d['gene']}\t{d['uniprot_id']}\t{d['pdb_file']}\n")

print("\nDone.")
