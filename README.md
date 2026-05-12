# LDT-TargetDB

**Ligand-Directed Transfer Covalent Radiopharmaceutical Target Database**

A structure-guided computational platform for systematic identification and prioritization of surface protein targets amenable to covalent radiopharmaceutical development.

## Features

- **2,886 human surface proteins** from SURFY in silico surfaceome
- **33 TCGA cancer types** × **30 GTEx normal tissues** expression profiling
- **AlphaFold-guided pocket detection** with nucleophile-specific scoring
- **LDT/NAS chemistry scoring** (LYS > TYR > HIS > SER)
- **Multi-strategy comparison**: LDT-NAS vs SuFEx-CTR vs traditional cysteine-targeted
- **Interactive web interface** with filtering, visualization, and data export

## Web Application

👉 **[LDT-TargetDB Live](https://ldt-targetdb.streamlit.app)** (deploying soon)

## Quick Start (Local)

```bash
pip install -r requirements.txt
streamlit run shiny_app/app.py
```

## Data Sources

| Layer | Source | Description |
|---|---|---|
| Surface proteome | SURFY (Bausch-Fluck et al., 2018 PNAS) | 2,886 predicted surface proteins |
| Expression | TCGA + GTEx (via UCSC Xena) | 33 cancer types, 30 normal tissues |
| Single-cell | TISCH2 | 190 scRNA-seq datasets |
| Structure | AlphaFold DB v6 | Protein 3D structures |
| Essentiality | DepMap 25Q4 | CRISPR Chronos scores |

## Citation

If you use LDT-TargetDB in your research, please cite:
*Manuscript in preparation.*

## License

MIT
