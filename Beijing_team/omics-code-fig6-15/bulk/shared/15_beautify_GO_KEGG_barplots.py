"""Redraw existing GO-BP and KEGG ORA results as publication bar plots.

No enrichment analysis is rerun. For each direction, the ten terms with the
smallest adjusted P values are displayed. Bar length is -log10(FDR), saturation
indicates whether FDR < 0.05, and labels report overlap gene counts.
"""

from __future__ import annotations

from pathlib import Path
import re
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\PCP\bulk\20260803\turn_20260802065502_636832")
FDR_THRESHOLD = 0.05
TOP_N = 10
BLUE = "#3B73B9"
RED = "#D64B40"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def clean_term(term: object) -> str:
    text = str(term).strip()
    text = re.sub(r"\s*\(GO:\d+\)\s*$", "", text)
    return "\n".join(textwrap.wrap(text, width=40, break_long_words=False))


def overlap_count(value: object) -> int:
    match = re.match(r"\s*(\d+)\s*/", str(value))
    return int(match.group(1)) if match else 0


def load_top_terms(path: Path, comparison: str, category: str, higher_in: str) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t")
    required = {"Term", "Overlap", "Adjusted P-value", "Combined Score", "Genes"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    data["FDR"] = pd.to_numeric(data["Adjusted P-value"], errors="coerce")
    data["Combined Score"] = pd.to_numeric(data["Combined Score"], errors="coerce")
    data = data.dropna(subset=["Term", "FDR"])
    data = data[data["FDR"] > 0].copy()
    data = data.sort_values(["FDR", "Combined Score"], ascending=[True, False]).head(TOP_N)
    if data.empty:
        raise ValueError(f"No usable enrichment terms in {path}")
    data["minus_log10_FDR"] = -np.log10(data["FDR"].clip(lower=1e-300))
    data["overlap_gene_count"] = data["Overlap"].map(overlap_count)
    data["comparison"] = comparison
    data["category"] = category
    data["higher_in"] = higher_in
    return data


def draw_bar_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    *,
    color: str,
    title: str,
    show_xlabel: bool,
) -> None:
    plot_data = data.sort_values(["FDR", "Combined Score"], ascending=[True, False]).copy()
    y = np.arange(len(plot_data))
    scores = plot_data["minus_log10_FDR"].to_numpy(dtype=float)
    alphas = np.where(plot_data["FDR"].to_numpy() < FDR_THRESHOLD, 0.92, 0.32)

    for position, score, alpha in zip(y, scores, alphas):
        ax.barh(
            position,
            score,
            height=0.68,
            color=color,
            alpha=float(alpha),
            edgecolor=color,
            linewidth=0.7,
        )

    threshold_x = -np.log10(FDR_THRESHOLD)
    max_score = max(float(scores.max()), threshold_x)
    right_limit = max_score * 1.20 + 0.18
    ax.axvline(threshold_x, color="#6D6D6D", linestyle=(0, (3, 2)), linewidth=0.85)
    ax.text(
        threshold_x,
        1.012,
        "FDR = 0.05",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="bottom",
        color="#555555",
        fontsize=7,
    )
    for position, score, count in zip(y, scores, plot_data["overlap_gene_count"]):
        ax.text(score + right_limit * 0.012, position, f"n={int(count)}", va="center", ha="left", fontsize=7, color="#333333")

    ax.set_yticks(y, labels=[clean_term(term) for term in plot_data["Term"]])
    ax.invert_yaxis()
    ax.set_xlim(0, right_limit)
    ax.set_title(title, fontsize=11, fontweight="bold", color=color, pad=13)
    ax.set_xlabel(r"$-\log_{10}$(FDR)" if show_xlabel else "")
    ax.tick_params(axis="y", length=0, pad=4, labelsize=8)
    ax.tick_params(axis="x", length=3, width=0.7)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#555555")


def save(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            output_base.with_suffix(f".{extension}"),
            dpi=600 if extension == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)


def draw_two_direction_figure(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    left_group: str,
    right_group: str,
    category_label: str,
    output_base: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.2))
    draw_bar_panel(
        axes[0], left, color=BLUE, title=f"Higher in {left_group}", show_xlabel=True
    )
    draw_bar_panel(
        axes[1], right, color=RED, title=f"Higher in {right_group}", show_xlabel=True
    )
    fig.text(0.5, 0.985, category_label, ha="center", va="top", fontsize=13, fontweight="bold")
    fig.subplots_adjust(left=0.20, right=0.985, bottom=0.10, top=0.88, wspace=0.72)
    save(fig, output_base)


def draw_combined_figure(
    go_left: pd.DataFrame,
    go_right: pd.DataFrame,
    kegg_left: pd.DataFrame,
    kegg_right: pd.DataFrame,
    *,
    left_group: str,
    right_group: str,
    output_base: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 11.0))
    draw_bar_panel(axes[0, 0], go_left, color=BLUE, title=f"GO BP — higher in {left_group}", show_xlabel=True)
    draw_bar_panel(axes[0, 1], go_right, color=RED, title=f"GO BP — higher in {right_group}", show_xlabel=True)
    draw_bar_panel(axes[1, 0], kegg_left, color=BLUE, title=f"KEGG — higher in {left_group}", show_xlabel=True)
    draw_bar_panel(axes[1, 1], kegg_right, color=RED, title=f"KEGG — higher in {right_group}", show_xlabel=True)
    fig.subplots_adjust(left=0.20, right=0.985, bottom=0.07, top=0.96, wspace=0.72, hspace=0.42)
    save(fig, output_base)


def process_comparison(
    *,
    folder: Path,
    comparison: str,
    left_group: str,
    right_group: str,
    output_stem: str,
) -> None:
    # Existing "down" tables represent the reference/left group; "up" tables
    # represent the comparison/right group in both analyses.
    go_left = load_top_terms(folder / "enrichr_down_GO_BP.tsv", comparison, "GO_BP", left_group)
    go_right = load_top_terms(folder / "enrichr_up_GO_BP.tsv", comparison, "GO_BP", right_group)
    kegg_left = load_top_terms(folder / "enrichr_down_KEGG.tsv", comparison, "KEGG", left_group)
    kegg_right = load_top_terms(folder / "enrichr_up_KEGG.tsv", comparison, "KEGG", right_group)

    selected = pd.concat([go_left, go_right, kegg_left, kegg_right], ignore_index=True)
    selected.to_csv(
        folder / f"barplot_GO_KEGG_{output_stem}_top10_source.tsv",
        sep="\t",
        index=False,
    )

    draw_two_direction_figure(
        go_left,
        go_right,
        left_group=left_group,
        right_group=right_group,
        category_label="GO Biological Process",
        output_base=folder / f"barplot_GO_BP_{output_stem}_beautified",
    )
    draw_two_direction_figure(
        kegg_left,
        kegg_right,
        left_group=left_group,
        right_group=right_group,
        category_label="KEGG Pathway",
        output_base=folder / f"barplot_KEGG_{output_stem}_beautified",
    )
    draw_combined_figure(
        go_left,
        go_right,
        kegg_left,
        kegg_right,
        left_group=left_group,
        right_group=right_group,
        output_base=folder / f"barplot_GO_KEGG_{output_stem}_combined_beautified",
    )
    counts = selected.groupby(["category", "higher_in"])["FDR"].apply(lambda x: int((x < FDR_THRESHOLD).sum()))
    print(f"{comparison}: selected {len(selected)} terms")
    print(counts.to_string())


def main() -> None:
    process_comparison(
        folder=ROOT / "02_PDC_vs_PDO" / "GO_KEGG",
        comparison="PDO_vs_PDC",
        left_group="PDC",
        right_group="PDO",
        output_stem="PDO_vs_PDC",
    )
    process_comparison(
        folder=ROOT / "03_Ulixertinib_vs_Ctrl" / "GO_KEGG",
        comparison="Ulixertinib_vs_Ctrl",
        left_group="Ctrl",
        right_group="Ulixertinib",
        output_stem="Ulixertinib_vs_Ctrl",
    )


if __name__ == "__main__":
    main()
