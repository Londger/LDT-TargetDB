Dear Editor,

We are pleased to submit our manuscript entitled "LDT-TargetDB: A Structure-Guided Covalent Radiopharmaceutical Target Database" for consideration in Briefings in Bioinformatics.

**Significance and Innovation**

Systematic target identification remains a fundamental challenge in covalent drug discovery. While covalent chemistry has transformed small-molecule therapeutics (e.g., KRAS G12C inhibitors, BTK inhibitors), no computational framework exists for systematically evaluating surface protein targets through the lens of covalent warhead chemistry. Existing tools (ImmunoTar, TCSA, DrugMap) are designed for immunotherapy targeting and lack structural chemistry dimensions essential for covalent probe design.

We present LDT-TargetDB, the first computational framework that integrates multi-omics data with 3D structural analysis and covalent chemistry scoring to systematically prioritize surface protein targets. Our key methodological innovations include:

1. **Integrated structural chemoproteomics pipeline.** We combine AlphaFold-guided binding pocket detection (pyKVFinder) with experimental pKa prediction (PROPKA) and solvent accessibility quantification (FreeSASA) to score individual nucleophilic residues (Lys>Cys>Tyr>Ser) within protein pockets — a capability absent from all existing surface-target tools.

2. **Covalent chemistry-aware scoring.** Unlike DrugMap (cysteine-only) or ImmunoTar (no structural chemistry), our framework explicitly models the ligand-directed transfer (LDT/NAS) reaction mechanism, incorporating residue type weights, pKa reactivity, SASA, and proximity-driven boost factors derived from published NASA chemistry literature.

3. **Multi-dimensional validation at scale.** We systematically validated our framework across 2,473 surface proteins using DepMap essentiality, Human Protein Atlas tissue expression, Open Targets disease associations, cross-species conservation, and TCGA mutation stability — demonstrating that top-ranked targets are non-essential, cancer-associated, highly conserved, and mutationally stable.

4. **Practical utility for experimentalists.** The framework provides an interactive web interface with real-time filtering, cancer-type-specific expression analysis (32 cancer types), and complete downloadable rankings, directly supporting experimental prioritization decisions.

In a 20-feature comparison, LDT-TargetDB implements 19/20 features, substantially exceeding DrugMap (11/20), ImmunoTar (9/20), and TCSA (8/20).

**Fit for Briefings in Bioinformatics**

Our manuscript addresses a recognized practical problem — how to computationally prioritize targets for covalent probe development — and provides a systematic solution backed by rigorous multi-dimensional validation. The framework is directly applicable to chemical biology, drug discovery, and structural bioinformatics. We believe this manuscript aligns with BIB's mission of providing indispensable methodological resources for experimental practitioners.

**Author Information**

Corresponding Author: Prof. Yiliang Li, Institute of Radiation Medicine, Chinese Academy of Medical Sciences & Peking Union Medical College, Tianjin 300192, China. Email: liyiliang@irm-cams.ac.cn.

This manuscript has not been published or submitted elsewhere. All authors have approved the submission.

Thank you for your consideration.

Sincerely,
Xinglong Zhou, Wenbin Hou, Yiliang Li
