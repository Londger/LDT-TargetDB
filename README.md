# LDT-TargetDB

**Ligand-Directed Transfer Covalent Radiopharmaceutical Target Database**

A structure-guided computational platform for systematic identification and prioritization of extracellular surface protein targets amenable to covalent radiopharmaceutical development.

## Features

- **2,799 unique human surface proteins** (SURFY in silico surfaceome after deduplication)
- **TCGA 32 cancer types** (10,535 samples; 9,186 primary tumors) × **GTEx 31 tissue sites** (7,862 samples) indication-aware expression profiling
- **AlphaFold-guided pocket detection** (pyKVFinder) with a refined extracellular-confined re-analysis of the top 200 targets
- **LDT/NASA chemistry scoring** with residue weights Lys (1.0) > Cys (0.3) > Tyr (0.2) > Ser (0.1)
- **Positive-unlabeled (PU) machine-learning prioritization**: deployed PU-v6 (44 features incl. tumor-microenvironment deconvolution and P2Rank ligandability), trained on a 54-target benchmark of known radioligand targets, all positives scored strictly out-of-fold — AUROC 0.830 (95% CI 0.785–0.868) with per-target bagged uncertainty
- **Indication-specific view** (per-cancer-type ITSI), annotated failure modes, clinical-evidence tiers and CT.gov trial counts
- **Interactive web interface** with filtering, 3D extracellular pocket visualization (py3Dmol), residue-level pKa/SASA tables, and CSV export

## Web Application

👉 **[LDT-TargetDB Live](https://ldt-targetdb.streamlit.app)**

## Quick Start (Local)

```bash
pip install -r requirements.txt
streamlit run shiny_app/app.py
```

## Data Sources

| Layer | Source | Description |
|---|---|---|
| Surface proteome | SURFY (Bausch-Fluck et al., 2018 PNAS) | 2,799 unique proteins after deduplication |
| Expression | TCGA + GTEx (via UCSC Xena) | 32 cancer types, 31 GTEx tissue sites |
| Structural | AlphaFold DB (v6) + pyKVFinder | 2,490 structures; refined EC-confined top-200 re-analysis |
| Protein evidence | Human Protein Atlas pathology IHC; UniProt N-glycosylation | Tumor protein levels; EC glycosylation site counts |
| Essentiality | DepMap 25Q4 | CRISPR Chronos scores |
| Ligandability | P2Rank 2.5.1 | Full-structure and extracellular-filtered runs |
| Clinical evidence | ClinicalTrials.gov (API v2), EU CT Register, literature | 54 benchmark targets; 867 radioligand trials |
| Single-cell | GSE162708, GSE72056 (targeted validation) | Blind-spot indications; systematic atlas integration planned |

## Citation

If you use LDT-TargetDB in your research, please cite:

*Zhou X, Hou W, Li Y. LDT-TargetDB: Systematic Prioritization of Extracellular Surface Targets for Covalent Radiopharmaceutical Development (under review).*

## License

MIT
