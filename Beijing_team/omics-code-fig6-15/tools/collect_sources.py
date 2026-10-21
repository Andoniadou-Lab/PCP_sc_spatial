"""Refresh this source-only snapshot from the adjacent PCP workspace."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from zipfile import ZipFile


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO.parent
BULK = Path("bulk/20260803/turn_20260802065502_636832")
SCRIPTS = BULK / "06_scripts"
FIG = Path("output/pdf/fig6_48pt_editable")
ZIP = FIG / "volcano_balanced_v2/Fig6_JM_FigS2_15_vector_PDFs_Arial48_balanced_v2.zip"

SOURCES = [
    (SCRIPTS / "01_sample_label_mapping.py", "bulk/pdc_pdo/01_sample_label_mapping.py"),
    (SCRIPTS / "02_quality_control.py", "bulk/pdc_pdo/02_quality_control.py"),
    (SCRIPTS / "03_PDC_vs_PDO_DE_figures.py", "bulk/pdc_pdo/03_PDC_vs_PDO_DE_figures.py"),
    (SCRIPTS / "04_pathway_gsea_GO_KEGG_PDC_PDO.py", "bulk/pdc_pdo/04_pathway_gsea_GO_KEGG_PDC_PDO.py"),
    (SCRIPTS / "05b_Ulixertinib_audit.py", "bulk/ulixertinib/05b_Ulixertinib_audit.py"),
    (SCRIPTS / "featureCounts_ulix.sh", "bulk/ulixertinib/featureCounts_ulix.sh"),
    (SCRIPTS / "05c_Ulixertinib_DE_figures.py", "bulk/ulixertinib/05c_Ulixertinib_DE_figures.py"),
    (SCRIPTS / "05f_Ulixertinib_GSEA_GO_KEGG.py", "bulk/ulixertinib/05f_Ulixertinib_GSEA_GO_KEGG.py"),
    (SCRIPTS / "06_gsea_full_lib.R", "bulk/shared/06_gsea_full_lib.R"),
    (SCRIPTS / "15_beautify_GO_KEGG_barplots.py", "bulk/shared/15_beautify_GO_KEGG_barplots.py"),
    (BULK / "_spreadsheet_work/gsea_selection/cache_gsea_rds.py", "bulk/shared/cache_gsea_rds.py"),
    (SCRIPTS / "19_plot_all_selected_gsea_terms.py", "bulk/shared/19_plot_all_selected_gsea_terms.py"),
    (FIG / "draw_go_kegg_48pt.py", "figures/draw_go_kegg_48pt.py"),
    (FIG / "draw_gsea_panels.py", "figures/draw_gsea_panels.py"),
    (FIG / "volcano_balanced_v2/draw_volcano_balanced_v2.py", "figures/volcano_balanced_v2/draw_volcano_balanced_v2.py"),
    (FIG / "volcano_balanced_v2/package_balanced_v2.py", "figures/volcano_balanced_v2/package_balanced_v2.py"),
    (FIG / "volcano_balanced_v2/README.md", "figures/volcano_balanced_v2/README.md"),
]

INPUTS = [
    Path("output/pdf/fig6_v3/data/volcano_plot_data.csv"),
    BULK / "02_PDC_vs_PDO/volcano/volcano_PDO_vs_PDC.csv",
    BULK / "02_PDC_vs_PDO/GO_KEGG/barplot_GO_KEGG_PDO_vs_PDC_top10_source.tsv",
    BULK / "03_Ulixertinib_vs_Ctrl/GO_KEGG/barplot_GO_KEGG_Ulixertinib_vs_Ctrl_top10_source.tsv",
    BULK / "05_tables/GSEA_selected_terms_plot_manifest.tsv",
    *[Path("output/pdf/fig6_v3/data") / f"{name}_curve.csv" for name in (
        "Ulixertinib_vs_Ctrl_GO0019827",
        "Ulixertinib_vs_Ctrl_KEGG_MAPK_hsa04010",
        "PDC_vs_PDO_GO0097150",
        "PDC_vs_PDO_HALLMARK_E2F_G2M",
        "PDC_vs_PDO_KEGG_MAPK_hsa04010",
    )],
]


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_csv(name: str, header: tuple[str, ...], rows: list[tuple]) -> None:
    with (REPO / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    sources = []
    for source, dest in SOURCES:
        payload = (ROOT / source).read_bytes()
        target = REPO / dest
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        sources.append((source.as_posix(), dest, len(payload), digest(payload)))
    write_csv("SOURCE_MANIFEST.csv", ("original_path", "repository_path", "bytes", "sha256"), sources)

    inputs = []
    for path in INPUTS:
        payload = (ROOT / path).read_bytes()
        inputs.append((path.as_posix(), len(payload), digest(payload)))
    write_csv("REQUIRED_INPUTS.csv", ("path_in_PCP_workspace", "bytes", "sha256"), inputs)

    with ZipFile(ROOT / ZIP) as archive:
        panels = [(entry.filename, entry.file_size, digest(archive.read(entry)))
                  for entry in archive.infolist() if entry.filename.lower().endswith(".pdf")]
    if len(panels) != 15:
        raise ValueError(f"Expected 15 PDFs in target ZIP, found {len(panels)}")
    write_csv("TARGET_PANELS.csv", ("path_in_target_zip", "bytes", "sha256"), panels)
    print(f"Copied {len(sources)} source files; indexed {len(inputs)} inputs and {len(panels)} target PDFs")


if __name__ == "__main__":
    main()
