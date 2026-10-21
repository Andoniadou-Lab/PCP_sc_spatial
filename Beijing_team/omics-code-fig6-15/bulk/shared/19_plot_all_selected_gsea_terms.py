"""Plot all stemness/proliferation/apoptosis/MAPK GSEA terms from saved results.

Downstream only: uses saved GSEA RDS objects, saved prerank lists, and saved VST
matrices. It does not rerun differential expression or GSEA permutations.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import pickle
import re
import textwrap
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r"E:\PCP\bulk\20260803\turn_20260802065502_636832")
CACHE = ROOT / "_spreadsheet_work" / "gsea_selection" / "gsea_rds_cache.pkl"
HALLMARK_CACHE = ROOT / "_spreadsheet_work" / "gsea_selection" / "MSigDB_Hallmark_2020.gmt.tsv"
MANIFEST_TSV = ROOT / "05_tables" / "GSEA_selected_terms_plot_manifest.tsv"
MANIFEST_JSON = ROOT / "05_tables" / "GSEA_selected_terms_plot_manifest.json"

BLUE = "#3775BA"
RED = "#B64342"
GRAY = "#6B7280"
LIGHT_BLUE = "#D9E8F6"
LIGHT_RED = "#F5DDDA"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9.5,
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)

# Direct pathway-name matches. These deliberately retain every matching result,
# regardless of nominal P value or FDR.
CATEGORY_PATTERNS = {
    "Stemness": re.compile(
        r"(?i)(STEMNESS|STEM[_ -]?CELL|PLURIPOTEN|SELF[_ -]?RENEWAL|EMBRYONIC[_ -]?STEM|IPS(?:C|[_ -]?CELL))"
    ),
    "Proliferation": re.compile(
        r"(?i)(PROLIFERAT|CELL[_ -]?CYCLE|E2F[_ -]?TARGET|G2M|G2[_ -]?M|MITOTIC[_ -]?SPINDLE)"
    ),
    "Apoptosis": re.compile(r"(?i)(APOPTOS|PROGRAMMED[_ -]?CELL[_ -]?DEATH)"),
    "MAPK": re.compile(
        r"(?i)(MAPK|MAP[_ -]?KINASE|ERK(?:1|2|5)?(?:_|\b)|RAF(?:_|\b)|MEK(?:1|2)?(?:_|\b)|JNK(?:_|\b)|P38(?:_|\b)|KRAS[_ -]?SIGNAL)"
    ),
}


def categories_for(*parts: object) -> list[str]:
    text = " ".join(str(p or "") for p in parts)
    return [name for name, pattern in CATEGORY_PATTERNS.items() if pattern.search(text)]


def safe_name(text: str, max_length: int = 110) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")
    return text[:max_length] or "term"


def format_p(value: float | None) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    if value == 0 or value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def read_hallmark_sets() -> dict[str, set[str]]:
    if not HALLMARK_CACHE.exists():
        url = (
            "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text"
            "&libraryName=MSigDB_Hallmark_2020"
        )
        with urllib.request.urlopen(url, timeout=60) as response:
            HALLMARK_CACHE.write_bytes(response.read())

    def normalize(text: str) -> str:
        return re.sub(r"[^A-Z0-9]+", "", text.upper())

    library: dict[str, set[str]] = {}
    for line in HALLMARK_CACHE.read_text(encoding="utf-8").splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        name = fields[0]
        genes = {g.strip() for g in fields[2:] if g.strip()}
        library[normalize(name)] = genes
    return library


def source_symbol_ranking(comparison: str) -> pd.DataFrame:
    if comparison == "PDC_vs_PDO":
        data = pd.read_csv(
            ROOT / "02_PDC_vs_PDO" / "differential_expression" / "all_genes_PDO_vs_PDC.tsv",
            sep="\t",
        )[["gene_symbol", "log2FC"]].rename(columns={"log2FC": "score"})
    else:
        data = pd.read_csv(
            ROOT / "03_Ulixertinib_vs_Ctrl" / "GSEA" / "prerank_input.tsv",
            sep="\t",
        ).iloc[:, :2]
        data.columns = ["gene_symbol", "score"]
    data = data.dropna(subset=["gene_symbol", "score"]).drop_duplicates("gene_symbol")
    return data.sort_values("score", ascending=False).reset_index(drop=True)


def expression_matrix(comparison: str) -> pd.DataFrame:
    if comparison == "PDC_vs_PDO":
        expression = pd.read_csv(ROOT / "01_quality_control" / "vst.tsv", sep="\t", index_col=0)
        annotation = pd.read_csv(
            ROOT / "02_PDC_vs_PDO" / "differential_expression" / "all_genes_PDO_vs_PDC.tsv",
            sep="\t",
        )
        mapping = (
            annotation.dropna(subset=["gene_symbol"])
            .drop_duplicates("ensembl")
            .set_index("ensembl")["gene_symbol"]
            .to_dict()
        )
        expression = expression.copy()
        mapping_base = {str(key).split(".")[0]: value for key, value in mapping.items()}
        expression["gene_symbol"] = [
            mapping.get(str(x), mapping_base.get(str(x).split(".")[0])) for x in expression.index
        ]
    else:
        expression = pd.read_csv(
            ROOT / "03_Ulixertinib_vs_Ctrl" / "vst.tsv", sep="\t", index_col=0
        ).T
        annotation = pd.read_csv(
            ROOT / "03_Ulixertinib_vs_Ctrl" / "count_matrix.metadata.tsv", sep="\t"
        )
        mapping = (
            annotation.dropna(subset=["gene_symbol"])
            .drop_duplicates("ensembl")
            .set_index("ensembl")["gene_symbol"]
            .to_dict()
        )
        expression = expression.copy()
        mapping_base = {str(key).split(".")[0]: value for key, value in mapping.items()}
        expression["gene_symbol"] = [
            mapping.get(str(x), mapping_base.get(str(x).split(".")[0])) for x in expression.index
        ]
    return (
        expression.dropna(subset=["gene_symbol"])
        .groupby("gene_symbol", sort=False)
        .mean(numeric_only=True)
    )


def map_entrez_to_symbol(ranked_entrez, source: pd.DataFrame) -> dict[str, str]:
    """Recover the R symbol-to-Entrez map by exact/nearest saved rank score."""
    source = source.sort_values("score").reset_index(drop=True)
    scores = source["score"].to_numpy(dtype=float)
    symbols = source["gene_symbol"].astype(str).to_numpy()
    mapping: dict[str, str] = {}
    values = ranked_entrez.values.astype(float)
    entrez = ranked_entrez.coords[ranked_entrez.dims[0]].values.astype(str)
    positions = np.searchsorted(scores, values)
    for gene_id, value, pos in zip(entrez, values, positions):
        candidates = []
        if pos < len(scores):
            candidates.append(pos)
        if pos > 0:
            candidates.append(pos - 1)
        best = min(candidates, key=lambda index: abs(scores[index] - value))
        if abs(scores[best] - value) <= 1e-10:
            mapping[gene_id] = symbols[best]
    return mapping


def running_curve(ranked_ids: np.ndarray, ranked_scores: np.ndarray, members: set[str]) -> pd.DataFrame:
    in_set = np.array([str(gene) in members for gene in ranked_ids], dtype=bool)
    n_hits = int(in_set.sum())
    if n_hits == 0 or n_hits == len(ranked_ids):
        raise ValueError(f"Invalid GSEA membership: {n_hits}/{len(ranked_ids)}")
    weighted_hits = np.abs(ranked_scores) * in_set
    hit_total = float(weighted_hits.sum())
    if hit_total == 0:
        weighted_hits = in_set.astype(float)
        hit_total = float(weighted_hits.sum())
    increments = np.where(in_set, weighted_hits / hit_total, -1.0 / (len(ranked_ids) - n_hits))
    return pd.DataFrame(
        {
            "rank": np.arange(len(ranked_ids), dtype=int),
            "gene_id": ranked_ids.astype(str),
            "rank_metric": ranked_scores.astype(float),
            "running_ES": np.cumsum(increments),
            "in_set": in_set,
        }
    )


def curve_from_pdc_csv(path: Path) -> tuple[pd.DataFrame, set[str]]:
    curve = pd.read_csv(path)
    curve["in_set"] = curve["in_set"].astype(str).str.lower().eq("true")
    curve = curve.rename(columns={"gene_symbol": "gene_id"})
    genes = set(curve.loc[curve["in_set"], "gene_id"].astype(str))
    return curve[["rank", "gene_id", "rank_metric", "running_ES", "in_set"]], genes


def curve_from_drug_dir(path: Path) -> tuple[pd.DataFrame, set[str]]:
    ranked = pd.read_csv(path / "prerank_data.rnk", sep="\t", header=None, names=["gene_id", "rank_metric"])
    fields = (path / "gene_sets.gmt").read_text(encoding="utf-8").strip().split("\t")
    members = set(fields[2:])
    curve = running_curve(
        ranked["gene_id"].astype(str).to_numpy(),
        ranked["rank_metric"].to_numpy(dtype=float),
        members,
    )
    return curve, members


def selected_terms(cache) -> list[dict]:
    hallmark_sets = read_hallmark_sets()
    records: list[dict] = []
    custom_terms = {
        "Stemness (Wong 2008 + Bhatt 2012)": "Stemness",
        "Proliferation (Hallmark E2F + G2M)": "Proliferation",
        "Apoptosis (Hallmark Apoptosis core)": "Apoptosis",
        "MAPK/ERK cascade (KEGG hsa04010)": "MAPK",
        "MAPK/ERK cascade (KEGG hsa04010 + Reactome RAF/MAPK)": "MAPK",
        "Cell cycle (GO 0007049)": "Proliferation",
        "Cell cycle (GO 0007049 + Reactome Cell Cycle)": "Proliferation",
    }

    for comparison in ("PDC_vs_PDO", "Ulixertinib_vs_Ctrl"):
        source = source_symbol_ranking(comparison)
        r_objects = cache[comparison]
        ranked_entrez = r_objects["ego"].geneList
        entrez_ids = ranked_entrez.coords[ranked_entrez.dims[0]].values.astype(str)
        entrez_scores = ranked_entrez.values.astype(float)
        symbol_map = map_entrez_to_symbol(ranked_entrez, source)

        # Custom/focused sets.
        gsea_folder = ROOT / ("02_PDC_vs_PDO" if comparison == "PDC_vs_PDO" else "03_Ulixertinib_vs_Ctrl") / "GSEA"
        focused = pd.read_csv(gsea_folder / "gsea_all_pathways.tsv", sep="\t")
        for term, category in custom_terms.items():
            match = focused.loc[focused["Term"].eq(term)]
            if match.empty:
                continue
            row = match.iloc[0]
            if comparison == "PDC_vs_PDO":
                file_map = {
                    "Stemness (Wong 2008 + Bhatt 2012)": "gsea_Stemness_Wong_2008___Bhatt_2012.csv",
                    "Proliferation (Hallmark E2F + G2M)": "gsea_Proliferation_Hallmark_E2F___G2M.csv",
                    "Apoptosis (Hallmark Apoptosis core)": "gsea_Apoptosis_Hallmark_Apoptosis_core.csv",
                    "MAPK/ERK cascade (KEGG hsa04010)": "gsea_MAPK_ERK_cascade_KEGG_hsa04010.csv",
                    "Cell cycle (GO 0007049)": "gsea_Cell_cycle_GO_0007049.csv",
                }
                curve, members = curve_from_pdc_csv(gsea_folder / file_map[term])
            else:
                dir_map = {
                    "Stemness (Wong 2008 + Bhatt 2012)": "_prerank_Stemness__Wong_2008___Bhatt_2012_",
                    "Proliferation (Hallmark E2F + G2M)": "_prerank_Proliferation__Hallmark_E2F___G2M_",
                    "Apoptosis (Hallmark Apoptosis core)": "_prerank_Apoptosis__Hallmark_Apoptosis_core_",
                    "MAPK/ERK cascade (KEGG hsa04010 + Reactome RAF/MAPK)": "_prerank_MAPK_ERK_cascade__KEGG_hsa04010___Reactome_RAF_MAPK_",
                    "Cell cycle (GO 0007049 + Reactome Cell Cycle)": "_prerank_Cell_cycle__GO_0007049___Reactome_Cell_Cycle_",
                }
                curve, members = curve_from_drug_dir(gsea_folder / dir_map[term])
            records.append(
                {
                    "comparison": comparison,
                    "category": category,
                    "source_sheet": "GSEA_Focused_All",
                    "library": "Focused",
                    "term_id": f"FOCUSED_{safe_name(term, 70)}",
                    "term_name": term,
                    "ES": float(row["ES"]),
                    "NES": float(row["NES"]),
                    "pvalue": float(row["NOM p-val"]),
                    "padj": float(row["FDR q-val"]),
                    "curve": curve,
                    "members_symbols": set(map(str, members)),
                }
            )

        # GO and KEGG from exact saved GSEA objects.
        for library, object_key in (("GO:BP", "ego"), ("KEGG", "ekegg")):
            obj = r_objects[object_key]
            result = obj.result
            for _, row in result.iterrows():
                term_id = str(row.get("ID", ""))
                term_name = str(row.get("Description", term_id))
                cats = categories_for(library, term_id, term_name)
                if not cats:
                    continue
                members = set(map(str, obj.geneSets.get(term_id, [])))
                if not members:
                    continue
                curve = running_curve(entrez_ids, entrez_scores, members)
                member_symbols = {symbol_map[g] for g in members if g in symbol_map}
                records.append(
                    {
                        "comparison": comparison,
                        "category": "/".join(cats),
                        "source_sheet": "GSEA_All_Libraries",
                        "library": library,
                        "term_id": term_id,
                        "term_name": term_name,
                        "ES": float(row.get("enrichmentScore", np.nan)),
                        "NES": float(row.get("NES", np.nan)),
                        "pvalue": float(row.get("pvalue", np.nan)),
                        "padj": float(row.get("p.adjust", np.nan)),
                        "curve": curve,
                        "members_symbols": member_symbols,
                    }
                )

        # Hallmark from saved fgsea results; membership is the same Hallmark 2020 library.
        hallmark = r_objects["efg"]
        source_ids = source["gene_symbol"].astype(str).to_numpy()
        source_scores = source["score"].to_numpy(dtype=float)

        def norm(text: str) -> str:
            return re.sub(r"[^A-Z0-9]+", "", text.upper().replace("HALLMARK_", ""))

        for _, row in hallmark.iterrows():
            term_id = str(row["pathway"])
            cats = categories_for("Hallmark", term_id)
            if not cats:
                continue
            members = hallmark_sets.get(norm(term_id), set())
            if not members:
                raise KeyError(f"Hallmark membership not found for {term_id}")
            curve = running_curve(source_ids, source_scores, members)
            records.append(
                {
                    "comparison": comparison,
                    "category": "/".join(cats),
                    "source_sheet": "GSEA_All_Libraries",
                    "library": "Hallmark",
                    "term_id": term_id,
                    "term_name": term_id,
                    "ES": float(row.get("ES", np.nan)),
                    "NES": float(row.get("NES", np.nan)),
                    "pvalue": float(row.get("pval", np.nan)),
                    "padj": float(row.get("padj", np.nan)),
                    "curve": curve,
                    "members_symbols": set(map(str, members)),
                }
            )

    # Stable category/library/term order and exact de-duplication.
    category_order = {"Stemness": 0, "Proliferation": 1, "Apoptosis": 2, "MAPK": 3}
    library_order = {"Focused": 0, "Hallmark": 1, "KEGG": 2, "GO:BP": 3}
    unique = {}
    for record in records:
        key = (record["comparison"], record["library"], record["term_id"])
        unique[key] = record
    return sorted(
        unique.values(),
        key=lambda item: (
            0 if item["comparison"] == "PDC_vs_PDO" else 1,
            min(category_order.get(cat, 9) for cat in item["category"].split("/")),
            library_order.get(item["library"], 9),
            item["term_name"],
        ),
    )


def score_values(expression: pd.DataFrame, members: set[str], comparison: str):
    present = sorted(set(expression.index.astype(str)).intersection(members))
    if not present:
        return None, None, present
    scores = expression.loc[present].mean(axis=0)
    if comparison == "PDC_vs_PDO":
        left_names = [f"PDC{i}" for i in range(1, 5)]
        right_names = [f"PDO{i}" for i in range(1, 5)]
    else:
        left_names = [f"Ctrl{i}" for i in range(1, 5)]
        right_names = [f"Ulixertinib{i}" for i in range(1, 5)]
    left = scores.reindex(left_names).to_numpy(dtype=float)
    right = scores.reindex(right_names).to_numpy(dtype=float)
    return left, right, present


def draw_curve(record: dict, output: Path) -> None:
    curve = record["curve"]
    x = curve["rank"].to_numpy(dtype=float)
    running = curve["running_ES"].to_numpy(dtype=float)
    metric = curve["rank_metric"].to_numpy(dtype=float)
    hits = curve["in_set"].to_numpy(dtype=bool)
    color = RED if record["NES"] >= 0 else BLUE
    if record["comparison"] == "PDC_vs_PDO":
        positive_group, negative_group = "PDO", "PDC"
    else:
        positive_group, negative_group = "Ulixertinib", "Ctrl"

    fig = plt.figure(figsize=(6.3, 5.0))
    grid = GridSpec(3, 1, figure=fig, height_ratios=[3.15, 0.48, 1.35], hspace=0.04)
    es_ax = fig.add_subplot(grid[0])
    hit_ax = fig.add_subplot(grid[1], sharex=es_ax)
    metric_ax = fig.add_subplot(grid[2], sharex=es_ax)

    es_ax.plot(x, running, color=color, lw=1.9)
    es_ax.fill_between(x, 0, running, color=color, alpha=0.11)
    es_ax.axhline(0, color=GRAY, lw=0.7)
    peak = int(np.argmax(np.abs(running)))
    es_ax.axvline(x[peak], color=color, lw=0.75, ls=(0, (2, 2)), alpha=0.75)
    title = "\n".join(textwrap.wrap(record["term_name"].replace("_", " "), width=58))
    es_ax.set_title(title, fontsize=10.5, fontweight="bold", pad=5)
    es_ax.set_ylabel("Running ES")
    es_ax.text(
        0.985,
        0.965,
        f"NES = {record['NES']:.2f}\nP = {format_p(record['pvalue'])}\nFDR = {format_p(record['padj'])}",
        transform=es_ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 2.5},
    )
    es_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    es_ax.spines["bottom"].set_visible(False)

    hit_ax.vlines(x[hits], 0, 1, color="#1F2937", lw=0.48)
    hit_ax.set_ylim(0, 1)
    hit_ax.set_yticks([])
    hit_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in hit_ax.spines.values():
        spine.set_visible(False)

    metric_ax.plot(x, metric, color="#666666", lw=0.55)
    metric_ax.fill_between(x, 0, metric, where=metric >= 0, color=RED, alpha=0.32, lw=0)
    metric_ax.fill_between(x, 0, metric, where=metric < 0, color=BLUE, alpha=0.32, lw=0)
    metric_ax.axhline(0, color=GRAY, lw=0.6)
    metric_ax.set_ylabel("Rank metric")
    metric_ax.set_xlabel("Rank in ordered gene list")
    metric_ax.text(0.01, -0.31, f"{positive_group} high", transform=metric_ax.transAxes,
                   ha="left", va="top", color=RED, fontsize=8)
    metric_ax.text(0.99, -0.31, f"{negative_group} high", transform=metric_ax.transAxes,
                   ha="right", va="top", color=BLUE, fontsize=8)
    metric_ax.set_xlim(x.min(), x.max())
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.18, top=0.88)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def draw_score(record: dict, left: np.ndarray, right: np.ndarray, n_genes: int, output: Path) -> float:
    if record["comparison"] == "PDC_vs_PDO":
        left_label, right_label = "PDC", "PDO"
    else:
        left_label, right_label = "Ctrl", "Ulixertinib"
    valid = np.isfinite(left) & np.isfinite(right)
    left, right = left[valid], right[valid]
    if len(left) >= 2:
        p_value = float(stats.ttest_rel(right, left).pvalue)
    else:
        p_value = np.nan

    fig, ax = plt.subplots(figsize=(5.1, 5.0))
    means = [float(np.mean(left)), float(np.mean(right))]
    sems = [float(stats.sem(left)) if len(left) > 1 else 0, float(stats.sem(right)) if len(right) > 1 else 0]
    ax.bar([0, 1], means, width=0.58, color=[LIGHT_BLUE, LIGHT_RED],
           edgecolor=[BLUE, RED], linewidth=1.2, zorder=1)
    for x_pos, mean, sem, color in zip([0, 1], means, sems, [BLUE, RED]):
        ax.errorbar(x_pos, mean, yerr=sem, fmt="none", ecolor=color,
                    elinewidth=1.3, capsize=4, capthick=1.1, zorder=2)
    jitters = np.linspace(-0.045, 0.045, len(left)) if len(left) else np.array([])
    for jitter, lv, rv in zip(jitters, left, right):
        ax.plot([jitter, 1 + jitter], [lv, rv], color="#8A8A8A", lw=1.15, alpha=0.74, zorder=3)
        ax.scatter(jitter, lv, s=42, color=BLUE, edgecolor="white", linewidth=0.65, zorder=4)
        ax.scatter(1 + jitter, rv, s=42, color=RED, edgecolor="white", linewidth=0.65, zorder=4)

    data_min = float(min(left.min(), right.min()))
    data_max = float(max(left.max(), right.max()))
    span = max(data_max - data_min, abs(data_max) * 0.04, 0.08)
    bracket_y = data_max + 0.12 * span
    bracket_h = 0.055 * span
    ax.plot([0, 0, 1, 1], [bracket_y, bracket_y + bracket_h, bracket_y + bracket_h, bracket_y],
            color="#333333", lw=0.9)
    ax.text(0.5, bracket_y + bracket_h + 0.025 * span, f"paired p = {format_p(p_value)}",
            ha="center", va="bottom", fontsize=9)
    title_text = textwrap.shorten(record["term_name"].replace("_", " "), width=86, placeholder="…")
    title = "\n".join(textwrap.wrap(title_text, width=43))
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=7)
    ax.set_xticks([0, 1], labels=[left_label, right_label])
    ax.set_ylabel("Gene-set expression score (mean VST)")
    ax.set_xlim(-0.48, 1.48)
    lower = 0.0 if data_min >= 0 else data_min - 0.22 * span
    ax.set_ylim(lower, bracket_y + bracket_h + 0.28 * span)
    ax.tick_params(axis="x", length=0)
    ax.spines["left"].set_color("#4B5563")
    ax.spines["bottom"].set_color("#4B5563")
    fig.subplots_adjust(left=0.19, right=0.985, bottom=0.13, top=0.82)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    return p_value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--scores-only", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    with CACHE.open("rb") as handle:
        cache = pickle.load(handle)
    records = selected_terms(cache)
    if args.limit is not None:
        records = records[: args.limit]

    counts = defaultdict(int)
    for record in records:
        counts[(record["comparison"], record["category"])] += 1
    print(f"Selected {len(records)} comparison-term records")
    for key, count in sorted(counts.items()):
        print(key, count)
    if args.inventory_only:
        return

    expression = {
        "PDC_vs_PDO": expression_matrix("PDC_vs_PDO"),
        "Ulixertinib_vs_Ctrl": expression_matrix("Ulixertinib_vs_Ctrl"),
    }
    manifest = []
    for index, record in enumerate(records, start=1):
        comparison_folder = "02_PDC_vs_PDO" if record["comparison"] == "PDC_vs_PDO" else "03_Ulixertinib_vs_Ctrl"
        out_dir = ROOT / comparison_folder / "GSEA" / "selected_all_related"
        stem = f"{record['library'].replace(':', '_')}_{safe_name(record['term_id'], 72)}"
        curve_path = out_dir / "curves" / f"{stem}_GSEA_curve.png"
        score_path = out_dir / "scores" / f"{stem}_score_dumbbell_bar.png"
        left, right, present = score_values(
            expression[record["comparison"]], record["members_symbols"], record["comparison"]
        )
        if left is None:
            print(f"SKIP no score genes: {record['comparison']} {record['library']} {record['term_id']}")
            continue
        if not args.scores_only:
            draw_curve(record, curve_path)
        paired_p = draw_score(record, left, right, len(present), score_path)
        manifest.append(
            {
                "plot_order": len(manifest) + 1,
                "comparison": record["comparison"],
                "category": record["category"],
                "source_sheet": record["source_sheet"],
                "library": record["library"],
                "term_id": record["term_id"],
                "term_name": record["term_name"],
                "ES": record["ES"],
                "NES": record["NES"],
                "pvalue": record["pvalue"],
                "padj": record["padj"],
                "n_genes_defined": len(record["members_symbols"]),
                "n_genes_scored": len(present),
                "paired_score_pvalue": paired_p,
                "curve_png": str(curve_path),
                "score_png": str(score_path),
            }
        )
        if index == 1 or index % 20 == 0 or index == len(records):
            print(f"PLOTTED {index}/{len(records)}", flush=True)

    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(MANIFEST_TSV, sep="\t", index=False)
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"MANIFEST {MANIFEST_TSV}")
    print(f"PLOTS {len(manifest)}")


if __name__ == "__main__":
    main()
