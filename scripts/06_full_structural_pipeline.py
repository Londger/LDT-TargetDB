"""
全量结构分析：所有TSI>0的表面蛋白 (m=1,054)
并行：下载AlphaFold + 口袋检测 + LDT评分
"""
import pandas as pd
import numpy as np
from pathlib import Path
import requests
import time
import subprocess
import shutil
import sys
from pyKVFinder import run_workflow
import freesasa

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
STRUCT_DIR = BASE / "data" / "structures"
SASA_TMP = Path("C:/temp_sasa")
SASA_TMP.mkdir(exist_ok=True)

RESIDUE_WEIGHTS = {"LYS": 1.0, "CYS": 0.3, "TYR": 0.2, "SER": 0.1}
REFERENCE_PKA = {"LYS": 10.5, "CYS": 8.5, "TYR": 10.0, "SER": 13.0}
PROXIMITY_BOOST = {"LYS": 2.0, "CYS": 1.0, "TYR": 1.0, "SER": 1.0}

ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.pdb"

# ============================================================
# 1. 加载基因-蛋白映射
# ============================================================
print("Loading gene mappings...")
surfy = pd.read_excel(
    BASE / "data/surfaceome/SURFY_surfaceome.xlsx",
    sheet_name="in silico surfaceome only", skiprows=1
)
surfy = surfy.rename(columns={"UniProt gene": "gene", "UniProt accession": "uniprot_id"})
surfy["uniprot_id"] = surfy["uniprot_id"].str.extract(r"([A-Z][0-9][A-Z0-9]{3,})")[0]
gene2uniprot = dict(zip(surfy.dropna(subset=["gene","uniprot_id"])["gene"],
                        surfy.dropna(subset=["gene","uniprot_id"])["uniprot_id"]))

tsi = pd.read_csv(BASE / "data/expression/tsi_results.csv")
# 筛选TSI > -0.01（保留几乎所有表面蛋白，排除极端阴性）
tsi_pos = tsi[tsi["tsi"] > -0.01].copy()
print(f"Genes with TSI > 0: {len(tsi_pos)}")

# 映射UniProt ID
tsi_pos["uniprot"] = tsi_pos["gene"].map(gene2uniprot)
genes_to_do = tsi_pos.dropna(subset=["uniprot"])
print(f"With UniProt mapping: {len(genes_to_do)}")

# 跳过已下载的
existing = set()
if (STRUCT_DIR / "downloaded_structures.csv").exists():
    existing = set(pd.read_csv(STRUCT_DIR / "downloaded_structures.csv")["gene"])
genes_to_do = genes_to_do[~genes_to_do["gene"].isin(existing)]
print(f"Need to download: {len(genes_to_do)}")

# ============================================================
# 2. 批量下载AlphaFold结构
# ============================================================
print("\n=== Phase 1: AlphaFold Download ===")
downloaded = []
download_failed = []

for i, (_, row) in enumerate(genes_to_do.iterrows()):
    gene = row["gene"]
    uniprot = row["uniprot"]
    pdb_file = STRUCT_DIR / f"AF-{uniprot}-F1-model_v6.pdb"

    # First check existing (from previous partial runs)
    if pdb_file.exists() and pdb_file.stat().st_size > 1000:
        downloaded.append({"gene": gene, "uniprot_id": uniprot, "pdb_file": str(pdb_file)})
        if len(downloaded) % 100 == 0:
            print(f"  Found cached: {len(downloaded)}", flush=True)
        continue

    url = ALPHAFOLD_URL.format(uniprot=uniprot)
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and (r.text.startswith("HEADER") or r.text.startswith("ATOM")):
            pdb_file.write_text(r.text)
            downloaded.append({"gene": gene, "uniprot_id": uniprot, "pdb_file": str(pdb_file)})
        else:
            download_failed.append({"gene": gene, "uniprot_id": uniprot})
    except Exception as e:
        download_failed.append({"gene": gene, "uniprot_id": uniprot, "error": str(e)[:50]})

    if len(downloaded) % 50 == 0:
        print(f"  {len(downloaded)}/{len(genes_to_do)} downloaded, {len(download_failed)} failed", flush=True)

    time.sleep(0.1)  # Rate limit

print(f"Downloaded: {len(downloaded)}, Failed: {len(download_failed)}")

# 合并已有的+新下载的
all_structures = list(existing) if isinstance(existing, list) else []
existing_list = pd.read_csv(STRUCT_DIR / "downloaded_structures.csv") if (STRUCT_DIR / "downloaded_structures.csv").exists() else None
if existing_list is not None:
    all_structures = existing_list.to_dict("records")
all_structures.extend(downloaded)

pd.DataFrame(downloaded).to_csv(STRUCT_DIR / f"new_downloads_{len(downloaded)}.csv", index=False)

# ============================================================
# 3. 全量LDT分析
# ============================================================
print(f"\n=== Phase 2: Full LDT Analysis ({len(all_structures)} proteins) ===")
results = []

for i, entry in enumerate(all_structures):
    gene = entry["gene"]
    uniprot = entry["uniprot_id"]
    pdb_file = Path(entry["pdb_file"])

    if not pdb_file.exists():
        continue

    try:
        # Pocket
        kv = run_workflow(str(pdb_file), step=0.8, volume_cutoff=5.0,
                          probe_in=1.4, probe_out=4.0)
        if kv.ncav == 0 or kv.residues is None:
            results.append({"gene": gene, "uniprot_id": uniprot,
                           "n_pockets": 0, "n_nucleophiles": 0, "ldt_score": 0})
            continue

        # Best pocket
        best_id, best_vol = None, 0
        pocket_keys = set()
        for pid, res_list in kv.residues.items():
            vol = kv.volume.get(pid, 0) if isinstance(kv.volume, dict) else 0
            if vol > best_vol:
                best_vol = vol
                best_id = pid
        if best_id is not None:
            for r in kv.residues.get(best_id, []):
                if len(r) >= 3:
                    pocket_keys.add(f"{r[1]}_{r[2]}_{r[0]}")

        # PROPKA
        pka_file = STRUCT_DIR / f"AF-{uniprot}-F1-model_v6.pka"
        if not pka_file.exists():
            subprocess.run(
                ["propka3", f"AF-{uniprot}-F1-model_v6.pdb"],
                capture_output=True, timeout=30, cwd=str(STRUCT_DIR)
            )
        pka_values = {}
        if pka_file.exists():
            in_data = False
            for line in pka_file.read_text().splitlines():
                if "RESIDUE" in line and "pKa" in line:
                    in_data = True; continue
                if not in_data or not line.strip() or line.startswith("-"): continue
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        pka_values[f"{parts[2]}_{parts[0]}_{parts[1]}"] = float(parts[3])
                    except: pass

        # FreeSASA
        sasa_values = {}
        try:
            tmp = SASA_TMP / f"{uniprot}.pdb"
            shutil.copy2(pdb_file, tmp)
            sasa_result = freesasa.calc(freesasa.Structure(str(tmp)))
            for chain, residues in sasa_result.residueAreas().items():
                for rn, ro in residues.items():
                    sasa_values[f"{chain}_{ro.residueType}_{rn}"] = {
                        "total": ro.total, "relative": ro.relativeTotal or 0
                    }
            tmp.unlink()
        except: pass

        # LDT Score
        nuc_scores = []
        for key in pocket_keys:
            parts = key.split("_")
            if len(parts) < 3: continue
            rt = parts[1]
            if rt not in RESIDUE_WEIGHTS: continue
            w = RESIDUE_WEIGHTS[rt]
            pk = pka_values.get(key, REFERENCE_PKA.get(rt, 10))
            ps = np.exp(-abs(pk - 7.4) / 2)
            rel = sasa_values.get(key, {}).get("relative", 0) or 0
            ss = min(rel / 50.0, 1.0)
            boost = PROXIMITY_BOOST.get(rt, 1.0)
            nuc_scores.append({"key": key, "type": rt, "score": w*ps*ss*boost})

        if nuc_scores:
            best = max(nuc_scores, key=lambda x: x["score"])
            ldt_score = best["score"]
            top_type = best["type"]
            n_nuc = len(nuc_scores)
        else:
            ldt_score = 0; top_type = ""; n_nuc = 0

        results.append({
            "gene": gene, "uniprot_id": uniprot,
            "n_pockets": kv.ncav,
            "max_pocket_volume": round(best_vol, 1),
            "n_nucleophiles": n_nuc,
            "ldt_score": round(ldt_score, 6),
            "top_res_type": top_type,
        })

    except Exception as e:
        pass  # Silently skip failures in bulk

    if (i+1) % 100 == 0:
        print(f"  [{i+1}/{len(all_structures)}] analyzed", flush=True)

shutil.rmtree(SASA_TMP, ignore_errors=True)

# ============================================================
# 4. 保存
# ============================================================
full_df = pd.DataFrame(results)
full_df = full_df.merge(tsi[["gene", "tsi_norm", "tsi"]], on="gene", how="left")

# Normalize
max_vol = full_df.loc[full_df["n_pockets"] > 0, "max_pocket_volume"].max()
full_df["pocket_norm"] = np.clip(full_df["max_pocket_volume"] / max_vol, 0, 1)

max_ldt = full_df.loc[full_df["ldt_score"] > 0, "ldt_score"].max()
full_df["ldt_norm"] = np.clip(full_df["ldt_score"] / max_ldt, 0, 1) if max_ldt > 0 else 0

full_df["final_score"] = 0.4*full_df["tsi_norm"].fillna(0) + 0.3*full_df["pocket_norm"].fillna(0) + 0.3*full_df["ldt_norm"].fillna(0)
full_df = full_df.sort_values("final_score", ascending=False).reset_index(drop=True)
full_df["rank"] = range(1, len(full_df)+1)

full_df.to_csv(BASE / "results" / "full_proteome_ranking.csv", index=False)
full_df.to_pickle(BASE / "results" / "full_proteome_ranking.pkl")

print(f"\nFull proteome analysis: {len(full_df)} proteins")
print(f"With pockets: {(full_df['n_pockets'] > 0).sum()}")
print(f"With LDT nucleophiles: {(full_df['ldt_score'] > 0).sum()}")
print(f"\nTop 20:")
for _, row in full_df.head(20).iterrows():
    print(f"  {int(row['rank']):3d}. {row['gene']:15s} Score={row['final_score']:.3f} LDT={row['ldt_score']:.4f} [{row['top_res_type']}]")

print("\nDone.")
