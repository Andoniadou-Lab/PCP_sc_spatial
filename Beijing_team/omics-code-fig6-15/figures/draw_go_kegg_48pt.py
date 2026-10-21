"""Export Fig6 K and FigS2 B enrichment panels as editable vector PDFs.

The bars and numerical annotations come from the saved top-10 enrichment
tables. No enrichment analysis is rerun. All term labels are Arial at 48 pt.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(r"D:\projects\PCP")
BASE = ROOT / "bulk/20260803/turn_20260802065502_636832"
sys.path.insert(0, str(BASE / "_spreadsheet_work/pydeps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle


OUTPUT = Path(__file__).resolve().parent / "panels"
OUTPUT.mkdir(parents=True, exist_ok=True)

BLUE = "#3977B7"
RED = "#BD4B49"
INK = "#202630"
GREY = "#727983"
GRID = "#E4E7EB"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": 48,
        "pdf.fonttype": 42,
        "pdf.use14corefonts": False,
        "path.simplify": False,
        "savefig.facecolor": "white",
    }
)

RENDERER = RendererAgg(1, 1, 72)
TERM_FONT = FontProperties(family="Arial", size=48)
PAGE_WIDTH = 2500.0
TERM_RIGHT = 1410.0
TERM_MAX_WIDTH = 1270.0
PLOT_X = 1500.0
PLOT_WIDTH = 730.0
COUNT_X = 2330.0
TOP = 215.0
BOTTOM = 190.0


def text_width(s: str) -> float:
    return float(RENDERER.get_text_width_height_descent(s, TERM_FONT, False)[0])


def wrap_label(s: str) -> list[str]:
    words = s.split()
    lines: list[str] = []
    current = ""
    for word in words:
        proposal = f"{current} {word}" if current else word
        if text_width(proposal) <= TERM_MAX_WIDTH:
            current = proposal
        else:
            if current:
                lines.append(current)
            if text_width(word) > TERM_MAX_WIDTH:
                raise ValueError(f"Term has a word wider than the label column: {word}")
            current = word
    if current:
        lines.append(current)
    assert lines and all(text_width(line) <= TERM_MAX_WIDTH for line in lines)
    return lines


def load_rows(source: Path, category: str, group: str) -> list[dict]:
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t")
                if r["category"] == category and r["higher_in"] == group]
    rows.sort(key=lambda r: (float(r["FDR"]), -float(r["Combined Score"])))
    if len(rows) != 10:
        raise ValueError(f"Expected 10 terms in {source.name}/{category}/{group}, found {len(rows)}")
    for row in rows:
        row["term_clean"] = re.sub(r"\s*\(GO:\d+\)\s*$", "", row["Term"])
        row["lines"] = wrap_label(row["term_clean"])
        row["value"] = float(row["minus_log10_FDR"])
        row["fdr"] = float(row["FDR"])
        row["count"] = int(row["overlap_gene_count"])
        if not math.isclose(row["value"], -math.log10(row["fdr"]), rel_tol=1e-8):
            raise ValueError(f"FDR value mismatch: {row['Term']}")
        if int(row["Overlap"].split("/")[0]) != row["count"]:
            raise ValueError(f"Overlap count mismatch: {row['Term']}")
    return rows


def figure_text(fig, x: float, y: float, value: str, **kwargs):
    width, height = fig.get_size_inches() * 72
    return fig.text(x / width, y / height, value, **kwargs)


def draw_panel(figure: str, category: str, group: str, source: Path) -> Path:
    rows = load_rows(source, category, group)
    heights = [max(107.0, 68.0 * len(row["lines"]) + 42.0) for row in rows]
    page_height = TOP + sum(heights) + BOTTOM
    fig = plt.figure(figsize=(PAGE_WIDTH / 72, page_height / 72), dpi=72)
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.set_xlim(0, PAGE_WIDTH)
    canvas.set_ylim(0, page_height)
    canvas.axis("off")

    color = BLUE if group in ("Ctrl", "PDC") else RED
    pale = "#D5E2EF" if group in ("Ctrl", "PDC") else "#EAC9C8"
    chart_type = "GO BP" if category == "GO_BP" else "KEGG"
    title = f"{chart_type} - higher in {group}"
    figure_text(fig, 95, page_height - 60, title, family="Arial", size=57,
                weight="bold", color=color, va="top")
    figure_text(fig, 95, page_height - 137, "Top 10 terms by FDR",
                family="Arial", size=34, color=GREY, va="top")

    values = [r["value"] for r in rows]
    xmax = max(max(values) * 1.11, -math.log10(0.05) * 1.22)
    first_y = page_height - TOP
    axis_y = first_y - sum(heights) - 30
    top_grid = first_y + 30
    threshold_x = PLOT_X + -math.log10(0.05) / xmax * PLOT_WIDTH

    # Grid and FDR threshold remain vector paths in the PDF.
    max_tick = int(math.floor(xmax))
    tick_interval = 1 if xmax < 5 else 2
    for tick in range(0, max_tick + 1, tick_interval):
        xx = PLOT_X + tick / xmax * PLOT_WIDTH
        canvas.plot([xx, xx], [axis_y, top_grid], color=GRID, lw=1.25, zorder=0)
        figure_text(fig, xx, axis_y - 22, str(tick), family="Arial", size=31,
                    color=INK, ha="center", va="top")
    canvas.plot([threshold_x, threshold_x], [axis_y, top_grid],
                color=GREY, lw=2.0, ls=(0, (4, 4)), zorder=2)
    figure_text(fig, threshold_x, top_grid + 14, "FDR = 0.05",
                family="Arial", size=29, color=GREY, ha="center", va="bottom")

    y_top = first_y
    for row, row_h in zip(rows, heights):
        y = y_top - row_h / 2
        label = "\n".join(row["lines"])
        figure_text(fig, TERM_RIGHT, y, label, family="Arial", size=48,
                    color=INK, ha="right", va="center", linespacing=1.08)
        bar_width = row["value"] / xmax * PLOT_WIDTH
        canvas.add_patch(Rectangle((PLOT_X, y - 26), bar_width, 52,
                                   facecolor=color if row["fdr"] < 0.05 else pale,
                                   edgecolor="none", zorder=1))
        figure_text(fig, COUNT_X, y, str(row["count"]), family="Arial",
                    size=34, color=INK, ha="center", va="center")
        y_top -= row_h

    canvas.plot([PLOT_X, PLOT_X + PLOT_WIDTH], [axis_y, axis_y], color=GREY, lw=1.8)
    figure_text(fig, PLOT_X + PLOT_WIDTH / 2, axis_y - 105, "-log10(FDR)",
                family="Arial", size=35, color=INK, ha="center", va="center")
    figure_text(fig, COUNT_X, top_grid + 14, "n", family="Arial", size=34,
                color=INK, ha="center", va="bottom", style="italic")

    path = OUTPUT / f"{figure}_{category}_{group}.pdf"
    fig.savefig(path, dpi=72, facecolor="white", metadata={
        "Title": f"{figure} {chart_type} higher in {group}",
        "Creator": "Matplotlib; source-backed vector redraw",
    })
    plt.close(fig)
    print(f"{path.name}\t{PAGE_WIDTH:.0f}x{page_height:.0f}pt\t" +
          f"rows={len(rows)}\tmax_label_lines={max(len(r['lines']) for r in rows)}")
    return path


def main() -> None:
    sources = [
        ("Fig6_K", BASE / "03_Ulixertinib_vs_Ctrl/GO_KEGG/barplot_GO_KEGG_Ulixertinib_vs_Ctrl_top10_source.tsv",
         ("Ctrl", "Ulixertinib")),
        ("FigS2_B", BASE / "02_PDC_vs_PDO/GO_KEGG/barplot_GO_KEGG_PDO_vs_PDC_top10_source.tsv",
         ("PDC", "PDO")),
    ]
    for figure, source, groups in sources:
        for category in ("GO_BP", "KEGG"):
            for group in groups:
                draw_panel(figure, category, group, source)


if __name__ == "__main__":
    main()
