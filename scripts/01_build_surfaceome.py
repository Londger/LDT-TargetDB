"""
模块1：表面蛋白初始集合构建
来源：SURFY (in silico) + UniProt API (GO验证)
"""
import pandas as pd
import json
import time
import requests
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
DATA = BASE / "data" / "surfaceome"
DATA.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. 读取 SURFY 表面蛋白数据（已验证可用）
# ============================================================
print("Loading SURFY data...")
surfy = pd.read_excel(
    DATA / "SURFY_surfaceome.xlsx",
    sheet_name="in silico surfaceome only",
    skiprows=1
)

# 清理：提取基因名和UniProt ID
surfy = surfy.rename(columns={
    "UniProt gene": "gene",
    "UniProt accession": "uniprot_id",
    "UniProt description": "description",
}).dropna(subset=["gene"])

# 清理HYPERLINK格式
surfy["uniprot_id"] = surfy["uniprot_id"].str.extract(r"([A-Z][0-9][A-Z0-9]{3,})")[0]
surfy = surfy.drop_duplicates(subset="gene")
genes_surfy = set(surfy["gene"].dropna())

print(f"  SURFY surface proteins: {len(genes_surfy)}")

# ============================================================
# 2. UniProt API 获取质膜/表面蛋白注释
# ============================================================
def fetch_uniprot_surface_genes():
    """
    从UniProt检索人类质膜蛋白（GO:0005886）和细胞表面蛋白（GO:0009986）
    """
    surface_genes = set()
    go_terms = [
        ("GO:0005886", "plasma_membrane"),
        ("GO:0005887", "integral_plasma_membrane"),
        ("GO:0009986", "cell_surface"),
    ]

    for go_id, go_name in go_terms:
        print(f"  Fetching {go_name} ({go_id})...")
        url = "https://rest.uniprot.org/uniprotkb/search"
        cursor = None
        count = 0

        while True:
            params = {
                "query": f"go:{go_id}",
                "fields": "gene_primary",
                "format": "json",
                "size": 500,
            }
            # TSV格式用tab，JSON用json
            headers = {"Accept": "application/json"}
            if cursor:
                params["cursor"] = cursor

            resp = requests.get(url, params=params, headers=headers)
            if resp.status_code != 200:
                print(f"    HTTP {resp.status_code}")
                break

            data = resp.json()
            results = data.get("results", [])
            if not results:
                break

            for entry in results:
                genes = entry.get("genes", [])
                for g in genes:
                    name = g.get("geneName", {}).get("value", "")
                    if name:
                        surface_genes.add(name)
                        count += 1

            # Next page
            next_link = resp.links.get("next", {}).get("url")
            if not next_link:
                break
            # Extract cursor from next link
            if "cursor=" in next_link:
                cursor = next_link.split("cursor=")[-1].split("&")[0]
            else:
                break

            time.sleep(0.3)

        print(f"    Got {count} genes (total unique: {len(surface_genes)})")

    return surface_genes

# 尝试用缓存避免重复请求
cache_file = DATA / "uniprot_surface_genes_cache.json"
if cache_file.exists():
    print("Loading UniProt cache...")
    with open(cache_file) as f:
        go_surface_genes = set(json.load(f))
else:
    try:
        go_surface_genes = fetch_uniprot_surface_genes()
        with open(cache_file, "w") as f:
            json.dump(list(go_surface_genes), f)
    except Exception as e:
        print(f"  UniProt API failed ({e}), proceeding with SURFY-only")
        go_surface_genes = set()

print(f"  GO surface genes: {len(go_surface_genes)}")

# ============================================================
# 3. 构建置信层级
# ============================================================
records = []
for gene in genes_surfy | go_surface_genes:
    in_surfy = gene in genes_surfy
    in_go = gene in go_surface_genes

    if in_surfy and in_go:
        tier = "Tier1_SURFY_GO"
    elif in_surfy:
        tier = "Tier2_SURFY_only"
    elif in_go:
        tier = "Tier3_GO_only"
    else:
        continue

    records.append({
        "gene": gene,
        "confidence_tier": tier,
        "source_surfy": in_surfy,
        "source_go": in_go,
    })

df = pd.DataFrame(records)

# 合并SURFY详细信息
surfy_info = surfy[["gene", "uniprot_id", "description", "length"]].copy()
df = df.merge(surfy_info, on="gene", how="left")

# 统计
print("\n=== Surface Protein Collection ===")
for tier in ["Tier1_SURFY_GO", "Tier2_SURFY_only", "Tier3_GO_only"]:
    n = (df["confidence_tier"] == tier).sum()
    print(f"  {tier}: {n}")

# 保存
df.to_csv(DATA / "surface_protein_collection.csv", index=False)
df.to_pickle(DATA / "surface_protein_collection.pkl")

# 输出用于后续分析的核心基因列表（优先Tier1+Tier2）
core_genes = df[df["confidence_tier"].isin(["Tier1_SURFY_GO", "Tier2_SURFY_only"])]
core_genes["gene"].to_csv(DATA / "core_surface_genes.txt", index=False, header=False)

print(f"\nCore surface genes (Tier1+2): {len(core_genes)}")
print("Done.")
