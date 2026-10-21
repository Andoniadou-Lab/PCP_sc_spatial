"""Create editable, 48 pt Arial GSEA panels from the saved plot inputs.

The curve data and statistics were produced in the earlier Fig6 rebuild. This
script changes layout and typography only; it does not rerun GSEA.
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(r"D:\projects\PCP")
BASE = ROOT / "bulk/20260803/turn_20260802065502_636832"
SOURCE = ROOT / "output/pdf/fig6_v3/data"
OUT = Path(__file__).resolve().parent
PANELS = OUT / "panels"
PREVIEWS = OUT / "previews"
for directory in (PANELS, PREVIEWS, OUT / "mplconfig"):
    directory.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(OUT / "mplconfig")
sys.path.insert(0, str(BASE / "_spreadsheet_work/pydeps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter, MaxNLocator

BLUE = "#3977B7"
RED = "#BD4B49"
GREY = "#727983"
INK = "#202630"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": 35,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.edgecolor": "#626973",
        "axes.linewidth": 1.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "path.simplify": False,
        "savefig.facecolor": "white",
    }
)

SPECS = [
    ("Fig6_L", "Ulixertinib_vs_Ctrl", "GO:0019827", "Ulixertinib_vs_Ctrl_GO0019827", "Stem cell population maintenance", "L"),
    ("Fig6_M", "Ulixertinib_vs_Ctrl", "hsa04010", "Ulixertinib_vs_Ctrl_KEGG_MAPK_hsa04010", "MAPK signaling pathway", "M"),
    ("FigS2_C", "PDC_vs_PDO", "GO:0097150", "PDC_vs_PDO_GO0097150", "Neuronal stem cell population maintenance", "C"),
    ("FigS2_D", "PDC_vs_PDO", "FOCUSED_Proliferation_Hallmark_E2F_G2M", "PDC_vs_PDO_HALLMARK_E2F_G2M", "Proliferation (Hallmark E2F + G2M)", "D"),
    ("FigS2_E", "PDC_vs_PDO", "hsa04010", "PDC_vs_PDO_KEGG_MAPK_hsa04010", "MAPK signaling pathway", "E"),
]


def format_p(value: float) -> str:
    return "< 0.001" if value < 0.001 else f"= {value:.3f}"


def draw_one(spec: tuple, stats: pd.DataFrame) -> dict:
    filename, comparison, term_id, curve_stem, title, letter = spec
    row = stats.loc[stats.comparison.eq(comparison) & stats.term_id.eq(term_id)]
    assert len(row) == 1, (comparison, term_id)
    stat = row.iloc[0]
    curve = pd.read_csv(SOURCE / f"{curve_stem}_curve.csv")
    x = curve["rank"].to_numpy(dtype=float) + 1
    es = curve["running_ES"].to_numpy(dtype=float)
    metric = curve["rank_metric"].to_numpy(dtype=float)
    hits = curve["in_set"].astype(str).str.lower().eq("true").to_numpy()
    peak = int(np.argmax(np.abs(es)))
    assert math.isclose(abs(es[peak]), abs(float(stat.ES)), abs_tol=1e-7)
    color = RED if stat.NES >= 0 else BLUE

    # A spacious individual-page layout gives 48 pt labels enough room without
    # changing the data geometry or cutting off the left axis labels.
    fig = plt.figure(figsize=(24, 15), facecolor="white")
    fig.text(0.035, 0.965, letter, fontweight="bold", fontsize=58, va="top")
    fig.text(0.135, 0.962, title, fontweight="bold", fontsize=48, va="top")
    fig.text(
        0.135,
        0.874,
        f"NES = {stat.NES:.2f}     P {format_p(stat.pvalue)}     FDR {format_p(stat.padj)}",
        fontsize=38,
        va="top",
    )

    left, width = 0.155, 0.800
    es_ax = fig.add_axes([left, 0.508, width, 0.297])
    hit_ax = fig.add_axes([left, 0.434, width, 0.050])
    rank_ax = fig.add_axes([left, 0.192, width, 0.218])

    es_ax.plot(x, es, color=color, lw=4.2)
    es_ax.fill_between(x, es, 0, color=color, alpha=0.10)
    es_ax.axhline(0, color=GREY, lw=1.3)
    es_ax.axvline(x[peak], color=color, lw=1.8, ls=(0, (3, 3)), alpha=0.8)
    es_ax.set_ylabel("Running ES", fontsize=48, labelpad=18)
    es_ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    es_ax.tick_params(axis="y", labelsize=33, length=6, pad=11)
    es_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    es_ax.spines["bottom"].set_visible(False)
    es_ax.margins(y=0.08)

    hit_ax.vlines(x[hits], 0, 1, color="#67707A", lw=1.0, alpha=0.85)
    hit_ax.set_ylim(0, 1)
    hit_ax.axis("off")

    rank_ax.fill_between(x, metric, 0, where=metric >= 0, color=RED, alpha=0.27)
    rank_ax.fill_between(x, metric, 0, where=metric < 0, color=BLUE, alpha=0.30)
    rank_ax.plot(x, metric, color="#818992", lw=1.3)
    rank_ax.axhline(0, color="#818992", lw=1.2)
    rank_ax.set_ylabel("Rank metric", fontsize=48, labelpad=18)
    rank_ax.set_xlabel("Rank in ordered gene list", fontsize=48, labelpad=18)
    rank_ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    rank_ax.tick_params(axis="both", labelsize=33, length=6, pad=11)
    tick_step = 3000 if len(curve) < 19000 else 4000
    rank_ax.set_xticks(np.arange(0, len(curve) + 1, tick_step))
    rank_ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))

    for ax in (es_ax, hit_ax, rank_ax):
        ax.set_xlim(0, len(curve))
    positive, negative = (
        ("PDO", "PDC") if comparison == "PDC_vs_PDO" else ("Ulixertinib", "Ctrl")
    )
    fig.text(left, 0.047, f"{positive} high", fontsize=38, color=RED)
    fig.text(left + width, 0.047, f"{negative} high", fontsize=38, color=BLUE, ha="right")

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas_bounds = fig.bbox
    outside = []
    for item in fig.findobj(matplotlib.text.Text):
        if not item.get_visible() or not item.get_text():
            continue
        bounds = item.get_window_extent(renderer)
        if (
            bounds.x0 < canvas_bounds.x0 - 1
            or bounds.y0 < canvas_bounds.y0 - 1
            or bounds.x1 > canvas_bounds.x1 + 1
            or bounds.y1 > canvas_bounds.y1 + 1
        ):
            outside.append(item.get_text())
    if outside:
        raise RuntimeError(f"{filename} has text outside page: {outside}")

    pdf_path = PANELS / f"{filename}.pdf"
    png_path = PREVIEWS / f"{filename}.png"
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=90)
    plt.close(fig)
    return {
        "pdf": str(pdf_path),
        "preview": str(png_path),
        "title": title,
        "letter": letter,
        "comparison": comparison,
        "term_id": term_id,
        "NES": float(stat.NES),
        "P": float(stat.pvalue),
        "FDR": float(stat.padj),
        "ranked_genes": len(curve),
        "hits": int(hits.sum()),
        "major_font_pt": 48,
        "minor_font_pt": [33, 38],
        "outside_page_text": outside,
    }


def main() -> None:
    stats = pd.read_csv(BASE / "05_tables/GSEA_selected_terms_plot_manifest.tsv", sep="\t")
    manifest = [draw_one(spec, stats) for spec in SPECS]
    (OUT / "gsea_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    for rec in manifest:
        print("SAVED", rec["pdf"], flush=True)


if __name__ == "__main__":
    main()
