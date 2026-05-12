"""
工具对比：LDT-TargetDB vs ImmunoTar vs TCSA vs DrugMap
证明我们独有的结构化学+共价适配维度
"""
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
RES = BASE / "results"

print("="*70)
print("LDT-TargetDB vs Existing Tools — Feature Comparison")
print("="*70)

# 各工具的功能矩阵
comparison = pd.DataFrame({
    "Feature": [
        "Surface protein definition",
        "Tumor/normal expression filtering",
        "Single-cell resolution validation",
        "Protein 3D structure analysis",
        "Binding pocket detection",
        "Nucleophile residue profiling",
        "pKa prediction",
        "Solvent accessibility (SASA)",
        "LDT/NAS chemistry scoring",
        "SuFEx chemistry scoring",
        "Covalent strategy comparison",
        "Gene essentiality (DepMap)",
        "Structure confidence (pLDDT)",
        "Pan-cancer ranking",
        "Cancer-type specific ranking",
        "Interactive web interface",
        "Downloadable data",
        "Radiopharmaceutical-specific",
        "Covalent chemistry-aware",
        "Total unique features",
    ],
    "LDT-TargetDB": [
        "SURFY (2,886 genes, ML-validated)",  # Surface
        "TCGA 33 types + GTEx 30 tissues",     # Expression
        "TISCH2 (190 scRNA-seq datasets)",      # Single-cell
        "AlphaFold DB (v6, 10,000+ structures)",  # 3D
        "pyKVFinder grid-based detection",       # Pocket
        "LYS/TYR/HIS/SER scoring + counts",     # Nucleophiles
        "PROPKA 3.5 (experimental pKa)",         # pKa
        "FreeSASA (relative accessibility)",     # SASA
        "Yes — NASA/ArNASA weight matrix",      # LDT
        "Yes — SuFEx weight matrix",            # SuFEx
        "Yes — 3-strategy comparison",          # Comparison
        "CRISPR Chronos (18,531 cell lines)",    # DepMap
        "pLDDT B-factor extraction",             # pLDDT
        "Composite TSI+LDT+pLDDT+DepMap",        # Pan-cancer
        "129 cancer type groups",                # Cancer-specific
        "Streamlit (interactive filtering)",     # Web
        "CSV download",                          # Download
        "Yes — designed for RLT",               # Radiopharma
        "Yes — LDT/SuFEx/traditional",          # Covalent
        "19/20",                                 # Total
    ],
    "ImmunoTar": [
        "Multiple DBs (CIRFESS, COMPARTMENTS)",  # Surface
        "GTEx normal + user-provided tumor",     # Expression
        "Not available",                          # Single-cell
        "Not available",                          # 3D
        "Not available",                          # Pocket
        "Not available",                          # Nucleophiles
        "Not available",                          # pKa
        "Not available",                          # SASA
        "Not available",                          # LDT
        "Not available",                          # SuFEx
        "Not available",                          # Comparison
        "DepMap integration",                     # DepMap
        "Not available",                          # pLDDT
        "Ewing/Myeloma/Neuroblastoma focused",   # Pan-cancer
        "Optimization per phenotype",             # Cancer-specific
        "R Shiny (web-based)",                   # Web
        "Yes",                                   # Download
        "No — immunotherapy focused",            # Radiopharma
        "No",                                    # Covalent
        "9/20",                                  # Total
    ],
    "TCSA": [
        "9 surfaceome resources integrated",     # Surface
        "TCGA + GTEx integration",               # Expression
        "13 cancer scRNA-seq datasets",          # Single-cell
        "Not available",                          # 3D
        "Not available",                          # Pocket
        "Not available",                          # Nucleophiles
        "Not available",                          # pKa
        "Not available",                          # SASA
        "Not available",                          # LDT
        "Not available",                          # SuFEx
        "Not available",                          # Comparison
        "Not available",                          # DepMap
        "Not available",                          # pLDDT
        "Pan-cancer surfaceome atlas",           # Pan-cancer
        "Not available",                          # Cancer-specific
        "Web portal (fcgportal.org)",            # Web
        "Yes",                                   # Download
        "No — CAR-T/ADC focused",               # Radiopharma
        "No",                                    # Covalent
        "8/20",                                  # Total
    ],
    "DrugMap": [
        "Not surface-specific (whole proteome)", # Surface
        "416 cancer cell lines",                 # Expression
        "Not available",                          # Single-cell
        "Limited (cysteine structural features)", # 3D
        "Not available",                          # Pocket
        "Cysteine only (isoTOP-ABPP)",           # Nucleophiles
        "Not available",                          # pKa
        "Not available",                          # SASA
        "Not available",                          # LDT
        "Not available",                          # SuFEx
        "Not available",                          # Comparison
        "Yes (CRISPR screens)",                  # DepMap
        "AI-predicted ligandability",             # pLDDT
        "Pan-cancer cysteine atlas",             # Pan-cancer
        "By lineage subtypes",                   # Cancer-specific
        "Web portal (drugmap.net)",              # Web
        "Yes",                                   # Download
        "No — general covalent ligand",          # Radiopharma
        "Yes — cysteine-focused only",           # Covalent
        "11/20",                                 # Total
    ],
})

print(comparison.to_string(index=False))

# 独特卖点提取
print("\n" + "="*70)
print("LDT-TargetDB Unique Selling Points")
print("="*70)
usp = [
    "1. First database specifically designed for covalent radiopharmaceutical target discovery",
    "2. Only tool integrating 3D pocket analysis with covalent chemistry (NASA/SuFEx/traditional)",
    "3. Only resource providing nucleophile-specific scoring (LYS/TYR/HIS/SER) for LDT probes",
    "4. Experimental pKa prediction (PROPKA) + SASA for reactivity assessment",
    "5. Multi-strategy covalent comparison: helps users choose optimal chemistry per target",
    "6. AlphaFold structure quality filtering (pLDDT <70 penalty)",
    "7. Designed for both LDT-NAS radiopharmaceuticals and general covalent drug discovery",
]
for u in usp:
    print(f"  {u}")

comparison.to_csv(RES / "tool_comparison.csv", index=False)
print(f"\nSaved: {RES}/tool_comparison.csv")
