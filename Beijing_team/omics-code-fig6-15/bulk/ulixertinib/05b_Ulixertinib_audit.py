"""Step 5c — Build Ulixertinib vs Ctrl sample sheet, cross-check vendor counts.

USAGE:
  python 06_scripts/05b_Ulixertinib_audit.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

TURN = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
A2 = TURN / "03_Ulixertinib_vs_Ctrl"
FC_DIR = A2 / "featureCounts"

# Build canonical sample sheet (same schema as bulk2/sample_sheet.csv)
# Vendor convention:  p3/p4/p5/p6 = Ctrl, p3u/p4u/p5u/p6u = Ulixertinib
# Each patient has one Ctrl + one Ulixertinib (paired design)
sample_sheet = pd.DataFrame([
    {"sample": "Ctrl1",         "patient_id": 3, "group": "Ctrl",        "drug": "DMSO",  "sequencing_batch": "data3", "bam": "p3.bam"},
    {"sample": "Ulixertinib1",  "patient_id": 3, "group": "Ulixertinib", "drug": "Ulixertinib", "sequencing_batch": "data3", "bam": "p3u.bam"},
    {"sample": "Ctrl2",         "patient_id": 4, "group": "Ctrl",        "drug": "DMSO",  "sequencing_batch": "data3", "bam": "p4.bam"},
    {"sample": "Ulixertinib2",  "patient_id": 4, "group": "Ulixertinib", "drug": "Ulixertinib", "sequencing_batch": "data3", "bam": "p4u.bam"},
    {"sample": "Ctrl3",         "patient_id": 5, "group": "Ctrl",        "drug": "DMSO",  "sequencing_batch": "data3", "bam": "p5.bam"},
    {"sample": "Ulixertinib3",  "patient_id": 5, "group": "Ulixertinib", "drug": "Ulixertinib", "sequencing_batch": "data3", "bam": "p5u.bam"},
    {"sample": "Ctrl4",         "patient_id": 6, "group": "Ctrl",        "drug": "DMSO",  "sequencing_batch": "data3", "bam": "p6.bam"},
    {"sample": "Ulixertinib4",  "patient_id": 6, "group": "Ulixertinib", "drug": "Ulixertinib", "sequencing_batch": "data3", "bam": "p6u.bam"},
])

sample_sheet.to_csv(A2 / "sample_sheet.csv", index=False)
sample_sheet.to_csv(A2 / "sample_sheet.tsv", sep="\t", index=False)
print("Sample sheet written:")
print(sample_sheet)

# Cross-check vendor counts (audit only — do not use for analysis)
vendor = pd.read_csv(A2 / "vendor_gene_count.xls", sep="\t")
print(f"\nVendor counts shape: {vendor.shape}")
print(f"Vendor columns: {list(vendor.columns[:10])}")
print(f"Samples in vendor: {[c for c in vendor.columns if c in ['p3','p3u','p4','p4u','p5','p5u','p6','p6u']]}")
# Check vendor vs our featureCounts once finished
print("\nNOTE: cross-check vs our featureCounts runs as a sanity comparison at end of pipeline.")