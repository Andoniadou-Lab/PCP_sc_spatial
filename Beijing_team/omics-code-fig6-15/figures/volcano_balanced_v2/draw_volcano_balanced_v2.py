"""Draw Fig6 J and FigS2 A with uniformly 48 pt Arial and the old layout ratio.

The physical page is 3.55 times the original single-panel page. This keeps
48 pt labels visually similar to the earlier 13-14 pt labels when fitted to a
page. Source data, thresholds, groups and highlighted genes are unchanged.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(r"D:\projects\PCP")
BASE = ROOT / "bulk/20260803/turn_20260802065502_636832"
sys.path.insert(0, str(BASE / "_spreadsheet_work/pydeps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
SCALE = 3.55
FONT = 48
INK = "#202630"
BLUE = "#3977B7"
RED = "#BD4B49"
GREY = "#7C8490"
PALE = "#CAD3DA"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": FONT,
        "pdf.fonttype": 42,
        "path.simplify": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.9 * SCALE,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "savefig.facecolor": "white",
    }
)

J_LABELS = {
    "HMGA2": (-2.70, 21.7), "NT5E": (-2.67, 19.4),
    "FERMT1": (-3.15, 17.55), "KRT18": (-1.25, 16.0),
    "DUSP4": (-0.63, 14.65), "PLCXD2": (-3.04, 13.6),
    "RAPGEF3": (-3.07, 12.05), "LINC00342": (-0.52, 12.7),
    "CYTOR": (-0.58, 11.15), "SHC4": (-3.15, 9.85),
    "SPON2": (2.27, 21.7), "EPB41L4B": (1.78, 16.5),
    "GPRIN2": (0.8, 9.9), "CMKLR2": (2.5, 9.25),
    "SLC6A9": (0.3, 8.4), "PCLO": (0.25, 6.45),
    "FAT2": (0.7, 4.35), "GJA1": (2.65, 6.75),
    "ANGPTL7": (3.02, 11.1), "MMP13": (3.1, 4.1),
}


def base_figure(title: str, xlim: tuple[float, float], ylim: tuple[float, float],
                xlabel: str, *, s2: bool) -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(10.0 * SCALE, 8.0 * SCALE), facecolor="white")
    ax = fig.add_axes([0.11, 0.13, 0.84, 0.71])
    fig.text(0.08, 0.96, title, fontsize=FONT, fontweight="bold", va="top")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([-6, -4, -2, 0, 2, 4] if s2 else [-4, -2, 0, 2, 4])
    ax.set_yticks([0, 2, 4, 6, 8, 10, 12] if s2 else [0, 5, 10, 15, 20])
    ax.set_xlabel(xlabel, fontsize=FONT, labelpad=8 * SCALE)
    ax.set_ylabel("−log10(FDR)", fontsize=FONT, labelpad=8 * SCALE)
    ax.tick_params(axis="both", labelsize=FONT, length=4 * SCALE,
                   width=0.9 * SCALE, pad=3 * SCALE)
    ax.axhline(-math.log10(0.05), color=GREY,
               ls=(0, (3, 3)), lw=0.9 * SCALE, zorder=1)
    for boundary in (-1, 1):
        ax.axvline(boundary, color=GREY, ls=(0, (3, 3)),
                   lw=0.9 * SCALE, zorder=1)
    fig.text(0.11, 0.035, "Dashed lines: |log2FC| = 1; FDR = 0.05",
             fontsize=FONT, color=GREY)
    return fig, ax


def check_and_save(fig: plt.Figure, name: str, expected_genes: list[str]) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_items = [t for t in fig.findobj(match=Text) if t.get_visible() and t.get_text()]
    for t in text_items:
        if abs(t.get_fontsize() - FONT) > 0.01:
            raise AssertionError((name, t.get_text(), t.get_fontsize()))
        box = t.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.x1 + 1 or box.y1 > fig.bbox.y1 + 1:
            raise AssertionError((name, "off page", t.get_text(), box))
    labels_found = {t.get_text() for t in text_items}
    if not set(expected_genes).issubset(labels_found):
        raise AssertionError((name, "missing labels", set(expected_genes) - labels_found))
    output = HERE / name
    fig.savefig(output, format="pdf", metadata={
        "Title": name.removesuffix(".pdf"),
        "Subject": "Source-backed, fully vector, all visible text Arial 48 pt",
    })
    plt.close(fig)
    print(output)


def fig6_j() -> None:
    data = pd.read_csv(ROOT / "output/pdf/fig6_v3/data/volcano_plot_data.csv")
    data["y"] = -np.log10(data.padj.clip(lower=1e-300))
    fig, ax = base_figure("Ulixertinib vs Ctrl", (-4.7, 4.05), (-0.6, 23),
                          "log2 fold change (Ulixertinib / Ctrl)", s2=False)
    groups = [
        ("ns", PALE, 8, 0.55),
        ("up_in_Ctrl", BLUE, 18, 0.78),
        ("up_in_Ulixertinib", RED, 18, 0.82),
    ]
    for direction, color, area, alpha in groups:
        rows = data.loc[data.direction.eq(direction)]
        ax.scatter(rows.log2FoldChange, rows.y, s=area * SCALE**2,
                   color=color, alpha=alpha, edgecolors="none", zorder=2)
    handles = [Line2D([], [], marker="o", ls="", color=color,
                      markersize=7 * SCALE, label=label)
               for color, label in ((PALE, "Not significant"),
                                    (BLUE, "Ctrl higher"),
                                    (RED, "Ulixertinib higher"))]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=3, fontsize=FONT, frameon=False, handletextpad=0.3,
              columnspacing=0.7, borderaxespad=0)
    for gene, pos in J_LABELS.items():
        row = data.loc[data.gene_symbol.eq(gene)].sort_values("padj").iloc[0]
        px, py = float(row.log2FoldChange), float(row.y)
        ax.scatter([px], [py], s=36 * SCALE**2, facecolors="white",
                   edgecolors=INK, lw=0.8 * SCALE, zorder=4)
        ax.annotate(gene, xy=(px, py), xytext=pos, ha="center", va="center",
                    fontsize=FONT, fontweight="normal", annotation_clip=False,
                    zorder=6,
                    bbox={"facecolor": "white", "edgecolor": "none",
                          "pad": 0.6, "alpha": 0.97},
                    arrowprops={"arrowstyle": "-", "lw": 0.65 * SCALE,
                                "color": GREY, "shrinkA": 2 * SCALE,
                                "shrinkB": 3 * SCALE})
    check_and_save(fig, "Fig6_J_volcano_Arial48_balanced_v2.pdf", list(J_LABELS))


def figs2_a() -> None:
    data = pd.read_csv(BASE / "02_PDC_vs_PDO/volcano/volcano_PDO_vs_PDC.csv")
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna(subset=["gene_symbol", "log2FC", "FDR", "significance"]).copy()
    data["y"] = -np.log10(data.FDR.clip(lower=1e-300))
    fig, ax = base_figure("PDO vs PDC", (-6.6, 4.9), (-0.5, 13.3),
                          "log2 fold change (PDO / PDC)", s2=True)
    groups = [("ns", PALE, 9, 0.50),
              ("down_sig", BLUE, 18, 0.78),
              ("up_sig", RED, 18, 0.80)]
    for group, color, area, alpha in groups:
        rows = data.loc[data.significance.eq(group)]
        ax.scatter(rows.log2FC, rows.y, s=area * SCALE**2,
                   color=color, alpha=alpha, edgecolors="none", zorder=2)
    handles = [Line2D([], [], marker=marker, ls="", color=color,
                      markersize=7 * SCALE, label=label)
               for color, label, marker in (
                   (RED, "PDO higher", "o"),
                   (BLUE, "PDC higher", "o"),
                   (PALE, "Not significant", "o"),
                   (INK, "Target gene", "D"),
               )]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=4, fontsize=FONT, frameon=False, handletextpad=0.3,
              columnspacing=0.7, borderaxespad=0)
    targets = ["SOX2", "SOX9", "FLRT2", "STOX2", "MKI67",
               "TACSTD2", "CASP3", "BCL2", "MAPK3", "MAPK1"]
    rows = data.loc[data.gene_symbol.isin(targets)].drop_duplicates("gene_symbol")
    found = []
    for side, subset in (("left", rows.loc[rows.log2FC < 0]),
                         ("right", rows.loc[rows.log2FC >= 0])):
        subset = subset.sort_values(["y", "gene_symbol"])
        slots = np.linspace(0.8, 5.5, len(subset))
        for (_, row), slot_y in zip(subset.iterrows(), slots):
            gene = str(row.gene_symbol)
            found.append(gene)
            px, py = float(row.log2FC), float(row.y)
            ax.scatter([px], [py], s=42 * SCALE**2, color=INK, marker="D",
                       edgecolors="white", lw=0.6 * SCALE, zorder=5)
            ax.annotate(gene, xy=(px, py),
                        xytext=(-6.05 if side == "left" else 4.35, slot_y),
                        ha="left" if side == "left" else "right", va="center",
                        fontsize=FONT, fontweight="bold", annotation_clip=False,
                        bbox={"facecolor": "white", "edgecolor": "#D1D5DB",
                              "pad": 1.7, "alpha": 0.97, "lw": 0.5 * SCALE},
                        arrowprops={"arrowstyle": "-", "lw": 0.7 * SCALE,
                                    "color": GREY, "shrinkA": 3 * SCALE,
                                    "shrinkB": 3 * SCALE}, zorder=6)
    check_and_save(fig, "FigS2_A_volcano_Arial48_balanced_v2.pdf", found)


if __name__ == "__main__":
    fig6_j()
    figs2_a()
