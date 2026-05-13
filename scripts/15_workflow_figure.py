"""
Figure 0: Professional Workflow Diagram
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

FIG = Path("c:/Users/Longer/Documents/自写生信论文/results/figures")
FIG.mkdir(exist_ok=True)

fig, ax = plt.subplots(figsize=(18, 10))
ax.set_xlim(0, 18)
ax.set_ylim(0, 10)
ax.axis('off')

# Color palette
C_INPUT  = '#2c3e50'   # dark blue-gray for input
C_PROC   = '#e74c3c'   # red for processing
C_OUTPUT = '#27ae60'   # green for output
C_WEB    = '#8e44ad'   # purple for web
C_ARROW  = '#7f8c8d'   # gray arrows
C_BG     = '#f8f9fa'   # light BG
C_ACC1   = '#3498db'
C_ACC2   = '#e67e22'
C_ACC3   = '#1abc9c'
C_ACC4   = '#9b59b6'
C_ACC5   = '#f39c12'

ax.set_facecolor(C_BG)

def draw_stage(ax, x, y, w, h, title, items, color, title_color='white'):
    """Draw a pipeline stage box"""
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                           facecolor=color, edgecolor='none', alpha=0.95)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h - 0.25, title, ha='center', va='top',
            fontsize=9, fontweight='bold', color=title_color)
    for j, item in enumerate(items):
        ax.text(x + w/2, y + h - 0.55 - j*0.22, item, ha='center', va='top',
                fontsize=6.5, color=title_color, alpha=0.85)

def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=2))

# ===== Title =====
ax.text(9, 9.6, 'LDT-TargetDB Analysis Pipeline', ha='center', fontsize=20, fontweight='bold', color='#2c3e50')
ax.text(9, 9.15, 'Integrated multi-omics + structural chemoproteomics for covalent radiopharmaceutical target discovery',
        ha='center', fontsize=10, color='#7f8c8d')

# ===== Row 1: Five Input Layers =====
y1 = 7.2
box_w, box_h = 2.8, 1.6
spacing = 0.5
start_x = 0.8

layers = [
    (C_ACC1, "Surface Proteome", ["SURFY in silico surfaceome", "2,799 proteins", "93.5% accuracy"]),
    (C_ACC2, "Expression Profiling", ["TCGA 33 cancer types", "GTEx 30 normal tissues", "9,186 tumors | 7,862 normals"]),
    (C_ACC3, "Cell-Type Validation", ["Human Protein Atlas IHC", "1,199,675 entries", "Known marker database"]),
    (C_ACC4, "3D Structure", ["AlphaFold DB v6", "2,495 structures", "pLDDT quality control"]),
    (C_ACC5, "Covalent Chemistry", ["PROPKA pKa prediction", "FreeSASA accessibility", "Nucleophile profiling"]),
]

for i, (color, title, items) in enumerate(layers):
    x = start_x + i * (box_w + spacing)
    draw_stage(ax, x, y1, box_w, box_h, title, items, color)

# ===== Arrows down to Processing =====
y2 = 5.5
for i in range(5):
    x = start_x + i * (box_w + spacing) + box_w/2
    draw_arrow(ax, x, y1, x, y2 + 1.2)

# ===== Row 2: Processing Pipeline =====
steps = [
    ("TSI Calculation", "Tumor Specificity Index\nPancancer ranking"),
    ("Pocket Detection", "pyKVFinder grid-based\nBinding site identification"),
    ("LDT Scoring", "LYS/CYS/TYR/SER weights\npKa × SASA × Proximity"),
    ("Composite Ranking", "TSI + Pocket + LDT + pLDDT\n+ DepMap essentiality"),
    ("Multi-Dimension\nValidation", "HPA | Open Targets | Ensembl\nTCGA mutations | STRING PPI"),
]

for i, (title, desc) in enumerate(steps):
    x = start_x + i * (box_w + spacing)
    draw_stage(ax, x, y2, box_w, 1.0, title, [desc], C_PROC)
    if i < 4:
        ax.annotate('', xy=(x + box_w + 0.15, y2 + 0.5), xytext=(x + box_w, y2 + 0.5),
                    arrowprops=dict(arrowstyle='->', color=C_PROC, lw=2))

# ===== Arrow to Output =====
draw_arrow(ax, start_x + 4*(box_w+spacing) + box_w/2, y2, 9, 3.5)

# ===== Row 3: Outputs =====
y3 = 1.8
outputs = [
    (C_OUTPUT, "LDT-TargetDB", ["2,473 ranked targets", "LYS-dominant (87% TOP30)", "19/20 features vs competitors"]),
    (C_OUTPUT, "Cancer-Type Specific", ["32 cancer types mapped", "9,563 samples (91%)", "Per-cancer best targets"]),
    (C_OUTPUT, "Validation Suite", ["DepMap: 99% non-essential", "21/30 cancer-associated", "All high conservation"]),
]
out_w = 5.0
out_spacing = 0.6
out_start = 0.8
for i, (color, title, items) in enumerate(outputs):
    x = out_start + i * (out_w + out_spacing)
    draw_stage(ax, x, y3, out_w, 1.3, title, items, color)

# ===== Web Tool =====
ax.text(9, 0.5, 'ltd-targetdb.streamlit.app', ha='center', fontsize=11,
        fontweight='bold', color=C_WEB, style='italic')
ax.text(9, 0.2, 'Interactive filtering · Target detail lookup · Cancer-type expression · CSV download',
        ha='center', fontsize=8, color='#95a5a6')

# ===== Stats sidebar =====
stats = [
    "Stats",
    f"2,799 surface proteins",
    f"2,663 with RNA-seq data",
    f"2,490 AlphaFold structures",
    f"2,478 with pockets (99.5%)",
    f"2,378 LDT nucleophiles (95.5%)",
    f"2,473 enhanced ranking",
    f"32 cancer types | 1,208 cell lines",
    f"195 structures pLDDT-validated",
]
for i, s in enumerate(stats):
    fs = 9 if i == 0 else 7
    fw = 'bold' if i == 0 else 'normal'
    ax.text(17.5, 8.5 - i * 0.28, s, fontsize=fs, fontweight=fw, ha='right',
            color='#2c3e50' if i == 0 else '#7f8c8d', fontfamily='monospace')

plt.tight_layout(pad=0.5)
fig.savefig(FIG / "figure0_workflow.png", dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 0 saved at 600 DPI")
