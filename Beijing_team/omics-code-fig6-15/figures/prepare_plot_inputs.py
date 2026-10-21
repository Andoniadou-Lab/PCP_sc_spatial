"""Prepare only the volcano and five GSEA curve CSVs used by the target PDFs.

This focused extraction follows the calculations in the original
output/pdf/fig6_v3/rebuild_figures.py. It requires the saved
gsea_exact_inputs.json, which is an upstream analysis output.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


TERMS = (
    ("Ulixertinib_vs_Ctrl", "GO:0019827", "Ulixertinib_vs_Ctrl_GO0019827"),
    ("Ulixertinib_vs_Ctrl", "hsa04010", "Ulixertinib_vs_Ctrl_KEGG_MAPK_hsa04010"),
    ("PDC_vs_PDO", "GO:0097150", "PDC_vs_PDO_GO0097150"),
    ("PDC_vs_PDO", "FOCUSED_Proliferation_Hallmark_E2F_G2M", "PDC_vs_PDO_HALLMARK_E2F_G2M"),
    ("PDC_vs_PDO", "hsa04010", "PDC_vs_PDO_KEGG_MAPK_hsa04010"),
)


def running(ids: list[str], scores: list[float], members: list[str]) -> pd.DataFrame:
    hit = np.isin(np.asarray(ids, dtype=str), np.asarray(members, dtype=str))
    weights = np.abs(np.asarray(scores, dtype=float))
    hit_weight = weights[hit].sum()
    if not hit.any() or hit.all() or hit_weight == 0:
        raise ValueError("Invalid GSEA membership or scores")
    increment = np.where(hit, weights / hit_weight, -1 / (len(hit) - hit.sum()))
    return pd.DataFrame({
        "rank": np.arange(len(hit)),
        "gene_id": ids,
        "rank_metric": scores,
        "running_ES": np.cumsum(increment),
        "in_set": hit,
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pcp-root", type=Path, default=Path(r"D:\projects\PCP"))
    parser.add_argument("--output", type=Path, help="Defaults to PCP/output/pdf/fig6_v3/data")
    args = parser.parse_args()
    root = args.pcp_root.resolve()
    bulk = root / "bulk/20260803/turn_20260802065502_636832"
    out = args.output.resolve() if args.output else root / "output/pdf/fig6_v3/data"
    out.mkdir(parents=True, exist_ok=True)

    exact = json.loads((root / "output/pdf/fig6_v3/data/gsea_exact_inputs.json").read_text(encoding="utf-8"))
    stats = pd.read_csv(bulk / "05_tables/GSEA_selected_terms_plot_manifest.tsv", sep="\t")
    for comparison, term, stem in TERMS:
        if term.startswith("FOCUSED_"):
            curve = pd.read_csv(
                bulk / "02_PDC_vs_PDO/GSEA/gsea_Proliferation_Hallmark_E2F___G2M.csv"
            ).rename(columns={"gene_symbol": "gene_id"})
            curve["in_set"] = curve["in_set"].astype(str).str.lower().eq("true")
        else:
            item = exact[comparison][term]
            curve = running(item["ranked_ids"], item["scores"], item["members"])
        row = stats.loc[stats.comparison.eq(comparison) & stats.term_id.eq(term)]
        if len(row) != 1:
            raise ValueError(f"Expected one GSEA statistic for {comparison}/{term}")
        peak = float(curve.running_ES.iloc[np.abs(curve.running_ES.to_numpy()).argmax()])
        if not np.isclose(peak, float(row.iloc[0].ES), atol=1e-8):
            raise ValueError(f"GSEA running ES mismatch for {stem}")
        destination = out / f"{stem}_curve.csv"
        curve.to_csv(destination, index=False)
        print(destination)

    drug_de = bulk / "03_Ulixertinib_vs_Ctrl/differential_expression/all_genes_Ulixertinib_vs_Ctrl.tsv"
    volcano = pd.read_csv(drug_de, sep="\t").replace([np.inf, -np.inf], np.nan)
    volcano = volcano.dropna(subset=["log2FoldChange", "padj"]).copy()
    volcano["y"] = -np.log10(volcano.padj.clip(lower=1e-300))
    destination = out / "volcano_plot_data.csv"
    volcano.to_csv(destination, index=False)
    print(destination)


if __name__ == "__main__":
    main()
