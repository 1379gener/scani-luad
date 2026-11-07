#!/usr/bin/env python3
"""Fig 1 workflow schematic for SCI-1 (IJMS). Schematic only: no data values are drawn;
cohort sizes are given in the legend text, not in the figure."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

import os
OUT = os.environ.get("FIG_OUT", "figures")  # repository default: ./figures
import os
os.makedirs(OUT, exist_ok=True)

plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.set_xlim(0, 100); ax.set_ylim(0, 72); ax.axis("off")

INK, MUTED = "#222222", "#666666"
FILL_MAIN, FILL_OUT, FILL_SIDE = "#EEF3F8", "#F6F6F6", "#FFFFFF"
EDGE_MAIN, EDGE_OUT = "#4A6A8A", "#7A7A7A"


def box(x, y, w, h, title, body, fill, edge, dashed=False):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                       linewidth=1.1, edgecolor=edge, facecolor=fill,
                       linestyle=(0, (4, 3)) if dashed else "solid")
    ax.add_patch(p)
    ax.text(x + w / 2, y + h - 2.6, title, ha="center", va="top", fontsize=8.6,
            fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - 7.4, body, ha="center", va="top", fontsize=7.2,
            color=MUTED, linespacing=1.35)
    return (x, y, w, h)


def arrow(x1, y1, x2, y2, dashed=False):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10,
                                 linewidth=1.0, color="#555555",
                                 linestyle=(0, (4, 3)) if dashed else "solid"))


# Row 1 — single-cell reference -> signature
b1 = box(2, 52, 44, 17, "Single-cell RNA-seq (GSE131907)",
         "Fibroblasts as annotated by the original study\n(primary tumour, normal lung, metastases)",
         FILL_MAIN, EDGE_MAIN)
b2 = box(54, 52, 44, 17, "myCAF signature (50 genes)",
         "Myofibroblasts vs. other fibroblast states\n(FindAllMarkers; top 50 by log2 fold change)",
         FILL_MAIN, EDGE_MAIN)
arrow(46.5, 60.5, 53.5, 60.5)

# Row 2 — index
b3 = box(27, 30, 50, 15, "SCANI",
         "Z(myCAF ssGSEA) + Z(TIMP1) − Z(PDPN)\neach component standardized within cohort",
         FILL_MAIN, EDGE_MAIN)
arrow(76, 51.5, 62, 45.6)
ax.text(86, 47.5, "TIMP1 and PDPN\nselected in\nTCGA-LUAD", ha="center", va="center",
        fontsize=6.8, color=MUTED, linespacing=1.3)

# Row 3 — cohorts
b4 = box(2, 3, 30, 19, "TCGA-LUAD",
         "Discovery, stage II–IV\nCox (per SD), adjusted Cox,\ntime-dependent ROC, KM", FILL_OUT, EDGE_OUT)
b5 = box(35, 3, 30, 19, "GSE72094",
         "Independent validation\n(microarray platform)\nCox (per SD), adjusted Cox, KM", FILL_OUT, EDGE_OUT)
b6 = box(68, 3, 30, 19, "CPTAC-LUAD",
         "Stage independence\n(Spearman correlation\nwith pathological stage)", FILL_OUT, EDGE_OUT)
for xc in (17, 50, 83):
    arrow(50, 29.5, xc, 22.6)

# Side — spatial, descriptive only (not a SCANI component)
box(2, 30, 21, 15, "Visium (GSE307534)",
    "Descriptive stromal\nprogram; not a\nSCANI component", FILL_SIDE, EDGE_OUT, dashed=True)

plt.tight_layout(pad=0.3)
for ext, kw in (("tiff", dict(dpi=600, pil_kwargs={"compression": "tiff_lzw"})),
                ("png", dict(dpi=600)), ("pdf", {})):
    fig.savefig(os.path.join(OUT, f"Fig1.{ext}"), **kw)
print("Fig1 written to", OUT)
