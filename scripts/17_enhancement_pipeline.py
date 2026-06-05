"""
增强分析管线: 胞外可及性 + Ligandability + Benchmark + 正常组织风险
"""
import pandas as pd, numpy as np, json, requests, time
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
RES = BASE / "results"
RES.mkdir(exist_ok=True)

# ============================================================
# 1. Extracellular Accessibility Filter
# ============================================================
print("=== 1. Extracellular Accessibility ===")

# Load existing data
enhanced = pd.read_csv(RES / "enhanced_final_ranking.csv")
surfy = pd.read_excel(BASE / "data/surfaceome/SURFY_surfaceome.xlsx",
                       sheet_name="in silico surfaceome only", skiprows=1)
surfy = surfy.rename(columns={"UniProt gene": "gene", "UniProt accession": "uniprot_id"})
surfy["uniprot_id"] = surfy["uniprot_id"].str.extract(r"([A-Z][0-9][A-Z0-9]{3,})")[0]
gene2uniprot = dict(zip(surfy.dropna(subset=["gene","uniprot_id"])["gene"],
                        surfy.dropna(subset=["gene","uniprot_id"])["uniprot_id"]))

# Get UniProt topology for top 200 targets
top200 = enhanced.head(200)
topo_cache = RES / "uniprot_topology_cache.json"
if topo_cache.exists():
    with open(topo_cache) as f:
        topology = json.load(f)
else:
    topology = {}
    for i, gene in enumerate(top200["gene"]):
        uniprot = gene2uniprot.get(gene)
        if not uniprot: continue
        try:
            r = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{uniprot}",
                params={"fields": "ft_topo_dom,ft_transmem,ft_signal,ft_chain"},
                headers={"Accept": "application/json"},
                timeout=10
            )
            if r.status_code == 200:
                data = r.json()
                features = []
                for feat in data.get("features", []):
                    ftype = feat.get("type","")
                    if ftype in ["Topological domain","Transmembrane","Signal peptide","Signal","Chain"]:
                        features.append({
                            "type": ftype,
                            "description": feat.get("description",""),
                            "start": feat.get("location",{}).get("start",{}).get("value",0),
                            "end": feat.get("location",{}).get("end",{}).get("value",0),
                        })
                topology[gene] = features
        except: pass
        if (i+1) % 50 == 0: print(f"  {i+1}/200 queried")
        time.sleep(0.15)
    with open(topo_cache, "w") as f:
        json.dump(topology, f)

# Classify extracellular regions
def classify_region(pos, features):
    """Check if a sequence position is extracellular"""
    if not features: return "unknown"
    for feat in features:
        if feat["start"] <= pos <= feat["end"]:
            if feat["type"] == "Topological domain" and "Extracellular" in feat.get("description",""):
                return "extracellular"
            elif feat["type"] == "Topological domain" and "Cytoplasmic" in feat.get("description",""):
                return "cytoplasmic"
            elif feat["type"] == "Transmembrane":
                return "transmembrane"
            elif feat["type"] == "Signal peptide":
                return "signal_peptide"
    return "unknown"

# Calculate extracellular accessibility score per target
# EAS = fraction of protein that is extracellular (based on topology)
eas_results = []
for gene in top200["gene"]:
    feats = topology.get(gene, [])
    # Count extracellular vs cytoplasmic vs TM regions
    n_extracellular = sum(1 for f in feats if "Extracellular" in str(f.get("description","")))
    n_cytoplasmic = sum(1 for f in feats if "Cytoplasmic" in str(f.get("description","")))
    n_transmembrane = sum(1 for f in feats if f["type"] == "Transmembrane")

    # EAS score:
    # 1.0 = has extracellular domains (true surface protein)
    # 0.7 = no topology data but SURFY predicted (likely surface)
    # 0.3 = only cytoplasmic/TM domains found (unlikely for SURFY proteins)
    if n_extracellular > 0:
        eas = 1.0
    elif len(feats) == 0:
        eas = 0.7  # No UniProt topology, but SURFY says surface
    elif n_transmembrane > 0 and n_extracellular == 0:
        eas = 0.5  # Has TM but no annotated extracellular domain
    else:
        eas = 0.3

    eas_results.append({
        "gene": gene,
        "has_topology_data": len(feats) > 0,
        "n_extracellular_domains": n_extracellular,
        "n_transmembrane_domains": n_transmembrane,
        "eas": 1.0 if n_extracellular > 0 else (0.5 if len(feats) > 0 else 0.0),
    })

eas_df = pd.DataFrame(eas_results)
# Merge back
enhanced_upd = enhanced.merge(eas_df[["gene","has_topology_data","eas","n_extracellular_domains"]],
                               on="gene", how="left")
enhanced_upd["eas"] = enhanced_upd["eas"].fillna(0.5)
print(f"  Genes with extracellular domains: {(eas_df['eas']==1.0).sum()}/{len(eas_df)}")
print(f"  Mean EAS: {eas_df['eas'].mean():.2f}")

# ============================================================
# 2. Enhanced Scoring with Extracellular Filter
# ============================================================
print("\n=== 2. Enhanced Final Scoring ===")

# New composite: TSI + Pocket + LDT + pLDDT + DepMap + EAS + Normal Tissue Risk
# Add extracellular penalty: targets without extracellular domains get penalized
enhanced_upd["eas_penalty"] = np.where(enhanced_upd["eas"] < 0.5, 0.5, 1.0)

enhanced_upd["enhanced_score_v2"] = (
    0.25 * enhanced_upd["tsi_norm"] +
    0.20 * enhanced_upd["pocket_norm"] +
    0.20 * enhanced_upd["ldt_norm"] * enhanced_upd["eas_penalty"] +
    0.10 * (1 - enhanced_upd["plddt_below_70"].fillna(100) / 100) +
    0.05 * np.where(enhanced_upd["depmap_label"] == "Essential", 1,
             np.where(enhanced_upd["depmap_label"] == "Context-dependent", 0.7,
             np.where(enhanced_upd["depmap_label"] == "No data", 0, 0.3))) +
    0.20 * enhanced_upd["eas"]  # Extracellular accessibility bonus
)

enhanced_upd = enhanced_upd.sort_values("enhanced_score_v2", ascending=False).reset_index(drop=True)
enhanced_upd["rank_v2"] = range(1, len(enhanced_upd)+1)

print("  Top 20 with extracellular filter:")
for _, row in enhanced_upd.head(20).iterrows():
    print(f"  {row['rank_v2']:3d}. {row['gene']:15s} Score={row['enhanced_score_v2']:.3f} EAS={row['eas']:.1f}")

# ============================================================
# 3. Benchmark: Known Targets
# ============================================================
print("\n=== 3. Benchmark Against Known Targets ===")

known_targets = {"TACSTD2","MET","LY6E","FAP","FOLH1","ERBB2","EGFR","DLL3","CLDN6",
                 "SSTR2","MSLN","CD46","LRRC15","NECTIN4","STEAP1","CA9","CLDN18"}
known_set = {g for g in enhanced_upd["gene"] if g.upper() in known_targets}

# Enrichment analysis
def enrichment_at_top_k(df, known_set, score_col, k_pcts=[0.01,0.05,0.10]):
    results = {}
    N = len(df)
    for pct in k_pcts:
        k = max(1, int(N * pct))
        top_k = set(df.head(k)["gene"])
        hit = len(top_k & known_set)
        expected = len(known_set) * pct
        enrichment = hit / expected if expected > 0 else 0
        results[f"top_{int(pct*100)}%"] = {
            "hits": hit, "total_targets": len(known_set),
            "enrichment": round(enrichment, 2)
        }
    return results

# Compare models
tsi_ranked = enhanced_upd.sort_values("tsi_norm", ascending=False).reset_index(drop=True)
ldt_ranked = enhanced_upd.sort_values("ldt_norm", ascending=False).reset_index(drop=True)
score_v2_ranked = enhanced_upd.sort_values("enhanced_score_v2", ascending=False).reset_index(drop=True)

models = {
    "Expression-only (TSI)": tsi_ranked,
    "LDT Chemistry-only": ldt_ranked,
    "Full model (v2)": score_v2_ranked,
}

for model_name, ranked_df in models.items():
    enrich = enrichment_at_top_k(ranked_df, known_set, "enhanced_score_v2")
    print(f"\n  {model_name}:")
    for k, v in enrich.items():
        print(f"    {k}: {v['hits']}/{v['total_targets']} hits, enrichment={v['enrichment']}x")

# ============================================================
# 4. Normal Tissue Risk Score
# ============================================================
print("\n=== 4. Normal Tissue Risk Score ===")

# Critical organs for radiopharmaceutical toxicity
CRITICAL_ORGANS = {
    "kidney": 3.0, "liver": 2.0, "salivary gland": 2.0,
    "bone marrow": 2.0, "spleen": 1.5, "intestine": 1.5,
    "lung": 1.0, "heart muscle": 1.0, "brain": 0.5,
}

# Load HPA normal tissue data
hpa = pd.read_csv(BASE / "data/expression/normal_ihc_data.tsv.zip",
                   sep="\t", compression="zip")

# Calculate NTRS for top targets
ntrs_results = []
for _, row in enhanced_upd.head(200).iterrows():
    gene = row["gene"]
    gene_hpa = hpa[hpa["Gene name"].str.upper() == gene.upper()]

    risk_score = 0
    n_organs_high = 0
    for organ, weight in CRITICAL_ORGANS.items():
        organ_data = gene_hpa[gene_hpa["Tissue"].str.lower().str.contains(organ.lower(), na=False)]
        if len(organ_data) > 0:
            high_levels = (organ_data["Level"] == "High").sum()
            medium_levels = (organ_data["Level"] == "Medium").sum()
            organ_risk = weight * (high_levels * 1.0 + medium_levels * 0.5) / len(organ_data) if len(organ_data) > 0 else 0
            risk_score += organ_risk
            if high_levels > 0: n_organs_high += 1

    ntrs_results.append({
        "gene": gene,
        "ntrs": round(risk_score, 2),
        "n_organs_high": n_organs_high,
        "risk_level": "Low" if n_organs_high == 0 else ("Medium" if n_organs_high <= 2 else "High"),
    })

ntrs_df = pd.DataFrame(ntrs_results)
enhanced_upd = enhanced_upd.merge(ntrs_df, on="gene", how="left")

print(f"  Targets with High risk: {(ntrs_df['risk_level']=='High').sum()}")
print(f"  Targets with Low risk: {(ntrs_df['risk_level']=='Low').sum()}")

# ============================================================
# 5. Save Enhanced Dataset
# ============================================================
output_cols = ["gene","enhanced_score_v2","rank_v2","tsi_norm","ldt_norm","pocket_norm",
               "eas","ntrs","risk_level","top_res_type","n_nucleophiles",
               "n_extracellular_domains","depmap_label","mean_plddt"]
output_cols = [c for c in output_cols if c in enhanced_upd.columns]
enhanced_upd[output_cols].to_csv(RES / "enhanced_ranking_v2.csv", index=False)

print(f"\n=== Enhancement Complete ===")
print(f"Saved: enhanced_ranking_v2.csv ({len(enhanced_upd)} targets)")
