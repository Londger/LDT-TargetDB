"""
模块3+4：口袋检测 + LDT可转移性评分
pyKVFinder + PROPKA + FreeSASA
"""
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import sys
import shutil
import os
import tempfile
from pyKVFinder import run_workflow
import freesasa
import warnings
warnings.filterwarnings("ignore")

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
STRUCT_DIR = BASE / "data" / "structures"
RESULTS_DIR = BASE / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ============================================================
# 1. LDT/NAS化学参数
# ============================================================
# NASA/ArNASA残基反应性权重（基于Tamura 2018 Nat Commun; Thimaradka 2021 BMC; Kawano 2023 JACS）
# Lys: 首选靶标，酰胺键产物稳定（kL ~10⁴ M⁻¹s⁻¹）
# Cys: 可反应但硫酯键不稳定，易水解
# Tyr: 可反应但酚酯键不稳定
# Ser: 可反应但烷酯键不稳定
# His: 文献明确报告"no labeled product observed" — 不反应（移除）
RESIDUE_WEIGHTS = {"LYS": 1.0, "CYS": 0.3, "TYR": 0.2, "SER": 0.1}
REFERENCE_PKA = {"LYS": 10.5, "CYS": 8.5, "TYR": 10.0, "SER": 13.0}
# LDT proximity effect: Lys ε-NH₂在pH 7.4下99.9%质子化，邻位效应(EM boost ~10⁶)补偿
PROXIMITY_BOOST = {"LYS": 2.0, "CYS": 1.0, "TYR": 1.0, "SER": 1.0}

# ============================================================
# 2. 加载数据
# ============================================================
print("Loading data...")
struct_map = pd.read_csv(STRUCT_DIR / "downloaded_structures.csv")
tsi_df = pd.read_csv(BASE / "data" / "expression" / "tsi_results.csv")
tsi_lookup = dict(zip(tsi_df["gene"], tsi_df["tsi_norm"]))

surfy = pd.read_excel(
    BASE / "data" / "surfaceome" / "SURFY_surfaceome.xlsx",
    sheet_name="in silico surfaceome only", skiprows=1
)
surfy = surfy.rename(columns={"UniProt gene": "gene"})
gene_lookup = dict(zip(surfy["gene"], surfy["UniProt description"]))

print(f"Structures: {len(struct_map)}")

# FreeSASA Temp ASCII directory (avoid Chinese path)
SASA_TMP = Path("C:/temp_sasa")
SASA_TMP.mkdir(exist_ok=True)

# ============================================================
# 3. 批量分析
# ============================================================
results = []

for i, row in struct_map.iterrows():
    gene = row["gene"]
    uniprot = row["uniprot_id"]
    pdb_file = Path(row["pdb_file"])

    if not pdb_file.exists():
        continue

    if (i + 1) % 20 == 0:
        print(f"  [{i+1}/{len(struct_map)}] downloaded={len(results)}", flush=True)

    try:
        # --- 3a. pqcket detection ---
        kv = run_workflow(str(pdb_file), step=0.8, volume_cutoff=5.0,
                          probe_in=1.4, probe_out=4.0)

        if kv.ncav == 0 or kv.residues is None or len(kv.residues) == 0:
            results.append({"gene": gene, "uniprot_id": uniprot,
                           "n_pockets": 0, "n_nucleophiles": 0,
                           "ldt_score": 0, "top_res_type": "", "top_res_pka": np.nan})
            continue

        # --- 3b. find largest pocket ---
        best_id, best_vol = None, 0
        pocket_keys = set()
        for pid, res_list in kv.residues.items():
            vol = kv.volume.get(pid, 0) if isinstance(kv.volume, dict) else 0
            if vol > best_vol:
                best_vol = vol
                best_id = pid

        # Convert pyKVFinder residue format: ['42', 'A', 'LYS'] -> "A_LYS_42"
        if best_id is not None:
            for r in kv.residues.get(best_id, []):
                if len(r) >= 3:
                    pocket_keys.add(f"{r[1]}_{r[2]}_{r[0]}")

        # --- 3c. PROPKA pKa ---
        pka_file = STRUCT_DIR / f"AF-{uniprot}-F1-model_v6.pka"
        if not pka_file.exists():
            # Run from STRUCT_DIR so output goes there
            subprocess.run(
                ["propka3", f"AF-{uniprot}-F1-model_v6.pdb"],
                capture_output=True, timeout=30, cwd=str(STRUCT_DIR)
            )

        pka_values = {}
        if pka_file.exists():
            in_data = False
            for line in pka_file.read_text().splitlines():
                if "RESIDUE" in line and "pKa" in line:
                    in_data = True
                    continue
                if not in_data or not line.strip() or line.startswith("-"):
                    continue
                parts = line.strip().split()
                # Format: RES_TYPE RES_NUM CHAIN pKa ...
                if len(parts) >= 4:
                    try:
                        res_type = parts[0]   # e.g., "LYS"
                        res_num = parts[1]    # e.g., "42"
                        chain = parts[2]      # e.g., "A"
                        pka = float(parts[3])
                        pka_values[f"{chain}_{res_type}_{res_num}"] = pka
                    except (ValueError, IndexError):
                        continue

        # --- 3d. FreeSASA (use ASCII temp path) ---
        sasa_values = {}
        try:
            tmp_pdb = SASA_TMP / f"{uniprot}.pdb"
            shutil.copy2(pdb_file, tmp_pdb)
            structure = freesasa.Structure(str(tmp_pdb))
            sasa_result = freesasa.calc(structure)
            for chain, residues in sasa_result.residueAreas().items():
                for resnum, res_obj in residues.items():
                    key = f"{chain}_{res_obj.residueType}_{resnum}"
                    sasa_values[key] = {
                        "total": res_obj.total,
                        "relative": res_obj.relativeTotal or 0
                    }
            tmp_pdb.unlink()
        except Exception:
            pass

        # --- 3e. LDT nucleophile scoring ---
        nuc_scores = []
        for key in pocket_keys:
            parts = key.split("_")
            if len(parts) < 3:
                continue
            rt = parts[1]  # residue type (3-letter)
            if rt not in RESIDUE_WEIGHTS:
                continue

            w = RESIDUE_WEIGHTS[rt]
            pk = pka_values.get(key, REFERENCE_PKA.get(rt, 10))
            ps = np.exp(-abs(pk - 7.4) / 2)
            sasa = sasa_values.get(key, {}).get("relative", 0) or 0
            ss = min(sasa / 50.0, 1.0)
            boost = PROXIMITY_BOOST.get(rt, 1.0)
            nuc_scores.append({
                "key": key, "type": rt, "score": round(w * ps * ss * boost, 4),
                "pka": round(pk, 1), "sasa": round(sasa, 1)
            })

        if nuc_scores:
            best = max(nuc_scores, key=lambda x: x["score"])
            ldt_score = best["score"]
            top_type = best["type"]
            top_pka = best["pka"]
            top_sasa = best["sasa"]
            n_nuc = len(nuc_scores)
        else:
            ldt_score = 0
            top_type = ""
            top_pka = np.nan
            top_sasa = 0
            n_nuc = 0

        results.append({
            "gene": gene, "uniprot_id": uniprot,
            "n_pockets": kv.ncav,
            "max_pocket_volume": round(best_vol, 1),
            "n_nucleophiles": n_nuc,
            "ldt_score": ldt_score,
            "top_res_type": top_type,
            "top_res_pka": top_pka,
            "top_res_sasa": top_sasa,
        })

    except Exception as e:
        results.append({
            "gene": gene, "uniprot_id": uniprot,
            "n_pockets": -1, "n_nucleophiles": 0,
            "ldt_score": 0, "top_res_type": "", "top_res_pka": np.nan,
            "error": str(e)[:80]
        })

# Cleanup temp
shutil.rmtree(SASA_TMP, ignore_errors=True)

# ============================================================
# 4. 综合排序
# ============================================================
print("\n=== Final Ranking ===")

final = pd.DataFrame(results)

# Merge TSI
final["tsi_norm"] = final["gene"].map(tsi_lookup).fillna(0)

# Normalize pocket volume
vp = final[final["n_pockets"] > 0]
max_vol = vp["max_pocket_volume"].max() if len(vp) > 0 else 1
final["pocket_norm"] = np.clip(final["max_pocket_volume"] / max_vol, 0, 1)

# Normalize LDT score
vl = final[final["ldt_score"] > 0]
max_ldt = vl["ldt_score"].max() if len(vl) > 0 else 1
final["ldt_norm"] = np.clip(final["ldt_score"] / max_ldt, 0, 1)

# Final composite score
final["final_score"] = (
    0.4 * final["tsi_norm"] +
    0.3 * final["pocket_norm"] +
    0.3 * final["ldt_norm"]
)

# Add description
final["description"] = final["gene"].map(gene_lookup).fillna("")

final = final.sort_values("final_score", ascending=False).reset_index(drop=True)
final["rank"] = range(1, len(final) + 1)

# Print summary
print(f"Analyzed: {len(final)}")
print(f"With pockets: {(final['n_pockets'] > 0).sum()}")
print(f"With LDT nucleophiles: {(final['ldt_score'] > 0).sum()}")
print(f"\n   Top 30 LDT Radiopharmaceutical Targets:\n")
print(f"{'Rank':<5} {'Gene':<14} {'Score':<8} {'TSI':<7} {'Pocket':<7} {'LDT':<7} {'Nuc':<5} {'BestNuc'}")
print("-" * 75)
for _, row in final.head(30).iterrows():
    nuc = f"{row.get('top_res_type',''):4s} pKa={row.get('top_res_pka','')}" if pd.notna(row.get('top_res_pka')) else "-"
    print(f"{row['rank']:<5} {row['gene']:<14} {row['final_score']:<8.3f} "
          f"{row['tsi_norm']:<7.3f} {row['pocket_norm']:<7.3f} {row['ldt_norm']:<7.3f} "
          f"{row.get('n_nucleophiles',0):<5} {nuc}")

# Save
final.to_csv(RESULTS_DIR / "final_ranking.csv", index=False)
final.to_pickle(RESULTS_DIR / "final_ranking.pkl")
final.head(20).to_csv(RESULTS_DIR / "top20_targets.csv", index=False)

print(f"\nResults saved to: {RESULTS_DIR}")
print("Done.")
