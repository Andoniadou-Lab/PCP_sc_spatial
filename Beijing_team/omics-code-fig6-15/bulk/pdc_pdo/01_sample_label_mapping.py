"""Step 1 — Sample label audit & mapping.

Rule (project spec, 2026-08-02):
  PDC = single-cell-derived culture;  PDO = sphere / organoid culture.
  Original labels from the sequencing vendor were *flipped*:
      "single cell" → actually sphere/organoid  → relabel as PDO
      "sphere / organoid" → actually single cell → relabel as PDC
  Final canonical names: PDC1..PDC3 / PDO1..PDO3 (or as many real samples exist).

Outputs:
  00_metadata_audit/sample_label_mapping.csv
  00_metadata_audit/sample_label_mapping.xlsx
  00_metadata_audit/sample_label_mapping.pdf
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

TURN_DIR = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
META_DIR = TURN_DIR / "00_metadata_audit"
META_DIR.mkdir(parents=True, exist_ok=True)

# Sources of truth
UPSTREAM_META = Path("/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/count_matrix.metadata.csv")
UPSTREAM_RAW = Path("/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/count_matrix.raw.tsv")

src = pd.read_csv(UPSTREAM_META)
cols = pd.read_csv(UPSTREAM_RAW, nrows=0, sep="\t").columns.tolist()
sample_cols = [c for c in cols if c != "Geneid"]
assert set(sample_cols) == set(src["sample"].tolist()), "Sample columns mismatch"

# Apply the spec: single_cells label was *flipped* → biology is sphere/organoid → PDO
#                  spheres label     was *flipped* → biology is single-cell-derived → PDC
def map_bio(bio: str) -> str:
    bio = str(bio).strip().lower()
    if bio == "single_cells":
        return "PDO"
    if bio == "spheres":
        return "PDC"
    raise ValueError(f"Unknown biology {bio!r}")

src["biology_corrected"] = src["biology"].map(map_bio)

# Final canonical sample name within each biology, by patient order
src = src.sort_values(["biology_corrected", "patient", "sample"]).reset_index(drop=True)
canonical_idx = {b: 0 for b in src["biology_corrected"].unique()}
canonical_names = []
for _, r in src.iterrows():
    canonical_idx[r["biology_corrected"]] += 1
    canonical_names.append(f"{r['biology_corrected']}{canonical_idx[r['biology_corrected']]}")
src["final_sample_name"] = canonical_names

# Inclusion decisions: include every patient pair (we have 4 pairs, 8 samples) — both for PDC and PDO
src["include_in_analysis"] = "yes"
src["exclusion_reason"] = ""

# Map direction strings
src["comparison_PDO_vs_PDC"] = "PDO vs PDC"   # all comparisons use this direction

# Build final mapping table (the 8 columns required by the spec)
mapping = src.rename(columns={
    "sample": "original_sample_name",
    "group": "original_group",
    "biology": "original_biology",
    "biology_corrected": "corrected_biology",
    "final_sample_name": "final_sample_name",
    "patient": "patient_id",
    "batch": "sequencing_batch",
    "include_in_analysis": "include_in_analysis",
    "exclusion_reason": "exclusion_reason",
    "comparison_PDO_vs_PDC": "comparison_direction",
}).copy()

# Reorder to spec order
mapping = mapping[
    [
        "original_sample_name",
        "original_group",
        "original_biology",
        "corrected_biology",
        "final_sample_name",
        "patient_id",
        "sequencing_batch",
        "include_in_analysis",
        "exclusion_reason",
        "comparison_direction",
    ]
]

# Also save a re-naming dict (handy for downstream code)
rename_dict = dict(zip(mapping["original_sample_name"], mapping["final_sample_name"]))

out_csv = META_DIR / "sample_label_mapping.csv"
mapping.to_csv(out_csv, index=False)

# Excel
out_xlsx = META_DIR / "sample_label_mapping.xlsx"
with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
    mapping.to_excel(w, sheet_name="label_mapping", index=False)

# PDF (simple text table via reportlab)
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet

out_pdf = META_DIR / "sample_label_mapping.pdf"
styles = getSampleStyleSheet()
doc = SimpleDocTemplate(str(out_pdf), pagesize=landscape(A4))
elements = [Paragraph("PCP project — sample label mapping (rule: D=single-cells label was flipped → PDO; Q=spheres label was flipped → PDC)", styles["Title"])]
elements.append(Spacer(1, 6))
elements.append(Paragraph("Comparison direction for downstream differential expression: PDO vs PDC (PDO numerator).", styles["BodyText"]))
elements.append(Spacer(1, 6))

tdata = [list(mapping.columns)] + mapping.astype(str).values.tolist()
tbl = Table(tdata, repeatRows=1)
tbl.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8EEF7")),
    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
    ("FONTSIZE", (0, 0), (-1, -1), 7),
    ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ("ALIGN", (0, 0), (-1, -1), "LEFT"),
]))
elements.append(tbl)
doc.build(elements)

# Save rename dict
import json
with open(META_DIR / "rename_dict.json", "w") as f:
    json.dump(rename_dict, f, ensure_ascii=False, indent=2)

print("MAPPING WRITTEN:", out_csv)
print(mapping.to_string(index=False))
print("\nRename dict:", rename_dict)