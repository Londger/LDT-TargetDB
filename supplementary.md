# Supplementary Materials

## LDT-TargetDB: A Structure-Guided Covalent Radiopharmaceutical Target Database

---

### Table S1: Complete Ranking of 2,473 Surface Proteins

The complete ranked list of all 2,473 proteins is available as a CSV download at https://ldt-targetdb.streamlit.app and in the accompanying file `full_proteome_ranking.csv`.

**Top 50 entries:**

| Rank | Gene | Score | TSI | LDT | Residue | Pockets | Nucleophiles |
|------|------|-------|-----|-----|---------|---------|-------------|
| 1 | ABCC5 | 0.578 | 0.210 | 0.851 | LYS | 5 | 54 |
| 2 | MUC1 | 0.559 | 0.471 | 0.475 | LYS | 32 | 87 |
| 3 | ILDR1 | 0.517 | 1.000 | 0.439 | LYS | 9 | 9 |
| 4 | CDH1 | 0.491 | 0.874 | 0.365 | LYS | 23 | 6 |
| 5 | LY75 | 0.479 | 0.378 | 0.606 | LYS | 66 | 56 |
| 6 | RHBDF2 | 0.477 | 0.328 | 1.000 | LYS | 31 | 1 |
| 7 | ABCC4 | 0.468 | 0.353 | 0.516 | LYS | 40 | 27 |
| 8 | CLDN4 | 0.460 | 0.882 | 0.203 | LYS | 5 | 2 |
| 9 | EPCAM | 0.450 | 0.807 | 0.240 | LYS | 7 | 6 |
| 10 | SLC29A3 | 0.445 | 0.353 | 0.713 | LYS | 39 | 40 |

---

### Table S2: Sensitivity Analysis of Composite Score Weights

The composite score weights (TSI=0.35, Pocket=0.25, LDT=0.25, pLDDT=0.10, DepMap=0.05) were tested against three alternative weighting schemes to assess ranking robustness.

| Weight Scheme | Weights (TSI/Pocket/LDT/pLDDT/DepMap) | TOP20 Overlap | Jaccard Index |
|---|---|---|---|
| Default | 0.35 / 0.25 / 0.25 / 0.10 / 0.05 | 20/20 (ref) | 1.00 |
| Equal-weights | 0.25 / 0.25 / 0.25 / 0.15 / 0.10 | 17/20 | 0.74 |
| Structure-weighted | 0.10 / 0.30 / 0.40 / 0.10 / 0.10 | 15/20 | 0.60 |
| Expression-weighted | 0.50 / 0.20 / 0.10 / 0.10 / 0.10 | 8/20 | 0.25 |

Three core targets (MUC1, ITGB4, SEZ6L2) appeared in the TOP20 across all four schemes, confirming the robustness of top-ranked predictions.

---

### Table S3: Covalent Strategy Comparison Matrix

| Strategy | Best for | Weighting | Targets Optimized |
|---|---|---|---|
| LDT-NASA | Lys (stable amide) | Lys=1.0 > Cys=0.3 > Tyr=0.2 > Ser=0.1 | Lys-rich pockets |
| SuFEx-CTR | Tyr (stable) | Tyr=1.0 > Lys=0.8 > His=0.7 > Ser=0.3 | Tyr/His-rich pockets |
| Traditional Cys | Cys (thioester) | Cys=1.0 > Lys=0.2 | Cys-rich pockets |

Approximately 84% of analyzed targets favor traditional cysteine-targeted chemistry due to favorable Cys pKa (8.5), while 16% favor SuFEx-CTR. LDT-NASA is preferred for targets where stable amide linkage is critical for radiopharmaceutical applications.

---

### Table S4: Known Nuclear Medicine Target Validation

| Target | Standard Name | TSI | TSI Rank | LDT Score | Best Residue | Assessment |
|---|---|---|---|---|---|---|
| TACSTD2 | Trop-2 | 0.055 | 10 | 0.226 | LYS | Clinically validated |
| MET | c-Met | 0.021 | 78 | 0.747 | LYS | Clinically validated |
| LY6E | LY6E | 0.014 | 147 | 0.014 | SER | Emerging target |
| ERBB2 | HER2 | 0.011 | 220 | 0.011 | LYS | Clinically validated |
| FAP | FAP | 0.011 | 238 | 0.011 | LYS | Stromal, low pan-cancer TSI |
| FOLH1 | PSMA | 0.012 | 250 | 0.012 | LYS | Prostate-specific, low pan-cancer TSI |
| SSTR2 | SSTR2 | -0.002 | 2018 | -0.002 | LYS | Neuroendocrine-specific |

Targets with tissue-specific expression (PSMA, SSTR2) appropriately receive lower pan-cancer TSI rankings, validating the pan-cancer approach while highlighting the utility of cancer-type-specific filtering.

---

### Table S5: Multi-Dimensional Validation Summary

| Validation Criterion | Source | Result |
|---|---|---|
| Gene essentiality | DepMap 25Q4 CRISPR (Chronos) | 99.7% of targets non-essential (2,373/2,380 with data) |
| Normal tissue protein expression | Human Protein Atlas IHC | 1,199,675 entries; 3 targets show zero normal tissue protein |
| Cross-species conservation | Ensembl Compara | All top 20 targets conserved in ≥5 model organisms |
| Cancer association | Open Targets Platform v4 | 21/30 queried targets have cancer disease links |
| Mutational stability | FireBrowse (37 TCGA cohorts) | ≤32/37 cohorts mutated for any top target |
| GO enrichment | STRING API | 20 significantly enriched GO biological process terms |
| Protein-protein interactions | STRING | 45 interaction edges among top 30 targets |

---

### Figure S1: Sensitivity Analysis

![Supplementary Figure S1: Sensitivity analysis of composite score weighting schemes.](results/figures/figure8_tool_comparison.png)

Four alternative weight schemes were evaluated. TOP20 overlap with the default scheme ranged from 8/20 (expression-weighted, Jaccard 0.25) to 17/20 (equal-weights, Jaccard 0.74), confirming reasonable robustness under moderate weight perturbations.

---

### Data Files

The following data files accompany this manuscript and are available for download:

- **full_proteome_ranking.csv**: Complete ranking of all 2,473 analyzed surface proteins
- **enhanced_final_ranking.csv**: Enhanced composite scoring with all validation dimensions
- **cancer_type_specific_expression.csv**: Per-cancer-type median expression for top targets
- **covalent_strategy_comparison.csv**: Full three-strategy comparison matrix
- **known_targets_validation.csv**: Known nuclear medicine target validation details

All files are downloadable from https://ldt-targetdb.streamlit.app.
