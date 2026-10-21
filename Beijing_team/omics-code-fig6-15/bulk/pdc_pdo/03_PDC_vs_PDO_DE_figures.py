"""Step 3 — Analysis 1: PDO vs PDC differential expression + volcano + heatmaps + scores.

Comparison direction (fixed): PDO vs PDC (PDO numerator, PDC denominator).
Design: DESeq2 with paired design (~ patient + group), group factor = corrected_biology
(PDC = reference; PDO = numerator).

Outputs:
  02_PDC_vs_PDO/differential_expression/all_genes_PDO_vs_PDC.tsv
  02_PDC_vs_PDO/differential_expression/significant_genes_PDO_vs_PDC.tsv
  02_PDC_vs_PDO/differential_expression/target_genes_PDO_vs_PDC.tsv
  02_PDC_vs_PDO/volcano/volcano_PDO_vs_PDC.png/pdf/svg + .csv
  02_PDC_vs_PDO/heatmap/heatmap_full.png/pdf/svg + .csv
  02_PDC_vs_PDO/heatmap/heatmap_targeted.png/pdf/svg + .csv
  02_PDC_vs_PDO/gene_expression/<gene>.{png,pdf,svg,csv}    for each required gene
  02_PDC_vs_PDO/pathway_scores/scores_PDO_vs_PDC.tsv + figures
"""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA

TURN = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
QC = TURN / "01_quality_control"
A1 = TURN / "02_PDC_vs_PDO"
DE = A1 / "differential_expression"
VOL = A1 / "volcano"
HM = A1 / "heatmap"
GE = A1 / "gene_expression"
PS = A1 / "pathway_scores"
for p in [DE, VOL, HM, GE, PS]:
    p.mkdir(parents=True, exist_ok=True)

RENAME = json.load(open(TURN / "00_metadata_audit" / "rename_dict.json"))

meta = pd.read_csv(TURN / "00_metadata_audit" / "sample_metadata_with_final_names.tsv", sep="\t")
filt = pd.read_csv(QC / "filtered_counts.tsv", sep="\t", index_col=0).astype(int)
expr = pd.read_csv(QC / "vst.tsv", sep="\t", index_col=0)
gene_info = pd.read_csv("/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/gene_info.tsv", sep="\t")
gene_info["ensembl"] = gene_info["Geneid"].str.split(".").str[0]
gene_info = gene_info.drop_duplicates("ensembl").set_index("ensembl")
gene_map = gene_info["gene_name"].to_dict()

# ----- Map sample -> corrected biology and patient -----
meta["corrected_biology"] = meta["biology"].map({"single_cells": "PDO", "spheres": "PDC"})
# pydeseq2 expects samples in obs; build a clean metadata with sample as index
deseq_meta = meta.set_index("sample")[["patient", "corrected_biology"]].copy()
deseq_meta["patient"] = deseq_meta["patient"].astype(str)
deseq_meta["corrected_biology"] = pd.Categorical(deseq_meta["corrected_biology"], categories=["PDC", "PDO"])
# Reorder filt columns to match metadata
filt = filt[deseq_meta.index.tolist()]

# ----- DESeq2 via pydeseq2 -----
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(
    counts=filt.T,
    metadata=deseq_meta,
    design_factors=["patient", "corrected_biology"],
    refit_cooks=True,
    n_cpus=2,
)
dds.deseq2()

# Pull PDO vs PDC (numerator=PDO)
stat_res = DeseqStats(dds, contrast=["corrected_biology", "PDO", "PDC"], n_cpus=2)
stat_res.summary()
res_df = stat_res.results_df.copy()
res_df.index.name = "ensembl"
res_df = res_df.reset_index()
res_df["gene_symbol"] = res_df["ensembl"].map(gene_map).fillna(res_df["ensembl"])
# Rename fields for clarity
res_df = res_df.rename(columns={
    "baseMean": "baseMean",
    "log2FoldChange": "log2FC",
    "lfcSE": "SE",
    "stat": "Wald_statistic",
    "pvalue": "P_value",
    "padj": "FDR",
})
res_df["comparison"] = "PDO vs PDC"
res_df["direction"] = np.where(res_df["log2FC"] > 0, "up_in_PDO",
                               np.where(res_df["log2FC"] < 0, "down_in_PDO", "ns"))
# significance
res_df["significance"] = "ns"
res_df.loc[(res_df["FDR"] < 0.05) & (res_df["log2FC"] > 1), "significance"] = "up_sig"
res_df.loc[(res_df["FDR"] < 0.05) & (res_df["log2FC"] < -1), "significance"] = "down_sig"

# Order columns
res_df = res_df[["ensembl", "gene_symbol", "baseMean", "log2FC", "SE",
                 "Wald_statistic", "P_value", "FDR", "comparison",
                 "direction", "significance"]]
res_df = res_df.sort_values("P_value")
res_df.to_csv(DE / "all_genes_PDO_vs_PDC.tsv", sep="\t", index=False)
res_df.to_csv(DE / "all_genes_PDO_vs_PDC.csv", index=False)

sig = res_df[(res_df["FDR"] < 0.05) & (res_df["log2FC"].abs() >= 1)].copy()
sig.to_csv(DE / "significant_genes_PDO_vs_PDC.tsv", sep="\t", index=False)
sig.to_csv(DE / "significant_genes_PDO_vs_PDC.csv", index=False)

# ----- TARGET GENES -----
TARGETS = [
    ("Stemness", "SOX2"),
    ("Stemness", "SOX9"),
    ("Stemness", "FLRT2"),
    ("Stemness", "STOX2"),
    ("Proliferation", "MKI67"),
    ("Proliferation", "TROP2 (TACSTD2)"),
    ("Apoptosis", "CASP3"),
    ("Apoptosis", "BCL2"),
]
# Symbol used to lookup in counts (real HGNC symbol)
TARGET_LOOKUP_SYMBOLS = ["SOX2", "SOX9", "FLRT2", "STOX2", "MKI67", "TACSTD2", "CASP3", "BCL2"]
TARGET_GENE_SYMBOLS = TARGET_LOOKUP_SYMBOLS[:]
# gene_map symbol lookup
sym_to_ens = {}
for ens, sym in gene_map.items():
    sym_to_ens.setdefault(sym, []).append(ens)
target_rows = []
for (cat, label), lookup_sym in zip(TARGETS, TARGET_LOOKUP_SYMBOLS):
    ens = sym_to_ens.get(lookup_sym, [])
    if not ens:
        target_rows.append({
            "category": cat, "display_label": label, "gene_symbol": lookup_sym,
            "ensembl": "not_detected", "log2FC": np.nan, "SE": np.nan,
            "Wald_statistic": np.nan, "P_value": np.nan, "FDR": np.nan,
            "baseMean": np.nan, "detection_status": "not_detected",
            "low_expression_filtered": "yes" if lookup_sym not in sym_to_ens else "no",
            "duplicate_gene_symbol": "no" if len(sym_to_ens.get(lookup_sym, [])) <= 1 else f"yes ({len(sym_to_ens[lookup_sym])} ensembl ids)",
        })
    else:
        for e in ens:
            r = res_df[res_df["ensembl"] == e]
            if r.empty:
                target_rows.append({
                    "category": cat, "display_label": label, "gene_symbol": lookup_sym,
                    "ensembl": e, "log2FC": np.nan, "SE": np.nan,
                    "Wald_statistic": np.nan, "P_value": np.nan, "FDR": np.nan,
                    "baseMean": np.nan, "detection_status": "filtered_due_to_low_expression",
                    "low_expression_filtered": "yes", "duplicate_gene_symbol": "no",
                })
            else:
                rr = r.iloc[0]
                detected = "yes" if rr["baseMean"] > 0 else "filtered_due_to_low_expression"
                target_rows.append({
                    "category": cat, "display_label": label, "gene_symbol": lookup_sym,
                    "ensembl": e, "log2FC": rr["log2FC"], "SE": rr["SE"],
                    "Wald_statistic": rr["Wald_statistic"], "P_value": rr["P_value"],
                    "FDR": rr["FDR"], "baseMean": rr["baseMean"],
                    "detection_status": detected,
                    "low_expression_filtered": "no" if rr["baseMean"] > 0 else "yes",
                    "duplicate_gene_symbol": "no" if len(ens) <= 1 else f"yes ({len(ens)} ensembl ids)",
                })
tg = pd.DataFrame(target_rows)
tg.to_csv(DE / "target_genes_PDO_vs_PDC.tsv", sep="\t", index=False)
tg.to_csv(DE / "target_genes_PDO_vs_PDC.csv", index=False)
print("Target gene table:\n", tg.to_string(index=False))

# ----- VOLCANO -----
def volcano(df, base, fdr_thr=0.05, lfc_thr=1.0):
    df = df.copy()
    df["neg_log10_FDR"] = -np.log10(df["FDR"].fillna(1.0).clip(lower=1e-300))
    color_map = {"up_sig": "#E07A5F", "down_sig": "#3B82C4", "ns": "#B8B8B8"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for cls, sub in df.groupby("significance"):
        ax.scatter(sub["log2FC"], sub["neg_log10_FDR"], c=color_map[cls],
                   s=14, alpha=0.65, edgecolor="none", label=f"{cls} (n={len(sub)})")
    # thresholds
    ax.axhline(-np.log10(fdr_thr), color="black", lw=0.8, ls="--")
    ax.axvline(lfc_thr, color="black", lw=0.8, ls="--")
    ax.axvline(-lfc_thr, color="black", lw=0.8, ls="--")
    # label target genes
    label_df = df[df["gene_symbol"].isin(TARGET_GENE_SYMBOLS)].copy()
    # If a gene symbol has multiple Ensembl rows, keep the one with highest baseMean
    label_df = label_df.sort_values("baseMean", ascending=False).drop_duplicates("gene_symbol")
    # Adjust label positions to avoid overlap
    label_df = label_df.sort_values("log2FC").reset_index(drop=True)
    for i, r in label_df.iterrows():
        sym = r["gene_symbol"]
        # Find display label
        for cat, lab in TARGETS:
            if lab.split(" (")[0] == sym:
                disp = lab
                break
        else:
            disp = sym
        x = r["log2FC"]; y = r["neg_log10_FDR"]
        # manual nudge for crowded points
        nudges = {  # symbol -> (dx, dy)
            "SOX2": (0.2, 0.3), "SOX9": (-0.5, 0.4), "FLRT2": (-0.4, 0.4),
            "STOX2": (-0.3, 0.4), "MKI67": (0.3, 0.4), "TACSTD2": (0.3, 0.4),
            "CASP3": (-0.3, 0.4), "BCL2": (-0.3, 0.4),
        }
        dx, dy = nudges.get(sym, (0.3, 0.4))
        ax.annotate(disp, xy=(x, y), xytext=(x+dx, y+dy),
                    fontsize=10, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="grey"))
    ax.set_xlabel("log2 fold change (PDO vs PDC)")
    ax.set_ylabel("-log10(FDR)")
    n_up = (df["significance"] == "up_sig").sum()
    n_down = (df["significance"] == "down_sig").sum()
    ax.set_title(f"Volcano — PDO vs PDC\nup_in_PDO = {n_up}  |  up_in_PDC = {n_down}  "
                 f"(FDR<{fdr_thr}, |log2FC|≥{lfc_thr})")
    ax.legend(loc="upper left", title="class", frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(VOL / f"{base}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    # also un-labeled version
    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    for cls, sub in df.groupby("significance"):
        ax.scatter(sub["log2FC"], sub["neg_log10_FDR"], c=color_map[cls],
                   s=10, alpha=0.6, edgecolor="none")
    ax.axhline(-np.log10(fdr_thr), color="black", lw=0.8, ls="--")
    ax.axvline(lfc_thr, color="black", lw=0.8, ls="--")
    ax.axvline(-lfc_thr, color="black", lw=0.8, ls="--")
    ax.set_xlabel("log2 fold change (PDO vs PDC)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"Volcano — PDO vs PDC (unlabeled)")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(VOL / f"{base}_nolabel.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    # legend-only figure
    fig, ax = plt.subplots(figsize=(3.6, 1.2))
    handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=8,
                          label=f"{cls} (n={ (df['significance']==cls).sum() })")
               for cls, c in color_map.items()]
    ax.legend(handles=handles, loc="center", frameon=False, ncol=1)
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(VOL / f"{base}_legend.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    # save underlying data
    df[["gene_symbol", "ensembl", "log2FC", "FDR", "baseMean", "significance"]].to_csv(
        VOL / f"{base}.csv", index=False)
    return df

vdf = volcano(res_df, "volcano_PDO_vs_PDC")

# ----- HEATMAP -----
PALETTE = {"PDC": "#3B82C4", "PDO": "#E07A5F"}
patient_palette = {10: "#66C2A5", 12: "#FC8D62", 14: "#8DA0CB", 15: "#E78AC3"}
batch_palette = {"data1": "#A6CEE3", "data2": "#FDBF6F"}

# Build a per-sample metadata frame for heatmap annotations
hm_meta = meta.copy()
hm_meta["corrected_biology"] = hm_meta["biology"].map({"single_cells": "PDO", "spheres": "PDC"})
hm_meta["final"] = hm_meta["sample"].map(RENAME)
hm_meta = hm_meta.set_index("final").loc[expr.columns.tolist()]

# A) Full heatmap of all filtered genes (12,528) — top 1000 by variance for visibility
top_var = expr.var(axis=1).sort_values(ascending=False).head(1000).index
mat = expr.loc[top_var]
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)
# Column ordering: group then patient
col_order = sorted(expr.columns, key=lambda c: (hm_meta.loc[meta.set_index("final_sample_name").loc[c, "sample"], "corrected_biology"] if False else 0,
                                                  c))
# simpler ordering: PDC1, PDC2, PDC3, PDC4, PDO1, PDO2, PDO3, PDO4
col_order = ["PDC1","PDC2","PDC3","PDC4","PDO1","PDO2","PDO3","PDO4"]
col_order = [c for c in col_order if c in z.columns]
z = z[col_order]

# Annotations
col_ann = pd.DataFrame({
    "group": [ hm_meta.loc[c, "corrected_biology"] for c in z.columns ],
    "patient": [ int(hm_meta.loc[c, "patient"]) for c in z.columns ],
    "batch": [ hm_meta.loc[c, "batch"] for c in z.columns ],
}, index=z.columns)
col_colors_group = col_ann["group"].map(PALETTE)
col_colors_patient = col_ann["patient"].map(patient_palette)
col_colors_batch = col_ann["batch"].map(batch_palette)

g = sns.clustermap(z, col_cluster=False, row_cluster=True, cmap="RdBu_r", center=0,
                   figsize=(7.5, 11), xticklabels=True, yticklabels=False,
                   col_colors=[col_colors_group, col_colors_patient, col_colors_batch])
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("genes (top 1000 by variance)")
g.ax_heatmap.set_title("Heatmap — top 1000 variable genes (z-score, full heatmap)", pad=14)
# legend
from matplotlib.patches import Patch
handles = [Patch(facecolor=PALETTE["PDC"], label="PDC"), Patch(facecolor=PALETTE["PDO"], label="PDO")]
g.figure.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
g.savefig(HM / "heatmap_full.png", dpi=600, bbox_inches="tight")
g.savefig(HM / "heatmap_full.pdf", bbox_inches="tight")
g.savefig(HM / "heatmap_full.svg", bbox_inches="tight")
plt.close(g.figure)

# B) Publication main heatmap — significant DE genes
sig_ens = sig["ensembl"].tolist()
mat = expr.loc[expr.index.isin(sig_ens)][col_order]
# keep top 200 by significance for visibility if many
if len(mat) > 200:
    mat = mat.loc[sig.set_index("ensembl").loc[mat.index].sort_values("FDR").head(200).index]
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)
g = sns.clustermap(z, col_cluster=False, row_cluster=True, cmap="RdBu_r", center=0,
                   figsize=(7.5, 11), xticklabels=True, yticklabels=False,
                   col_colors=[col_colors_group, col_colors_patient, col_colors_batch])
g.ax_heatmap.set_title("Heatmap — top significant DE genes (FDR<0.05, |log2FC|≥1)", pad=14)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("significant DE genes")
g.figure.legend(handles=handles + [Patch(facecolor=patient_palette[p], label=f"patient {p}") for p in sorted(patient_palette)],
                loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
g.savefig(HM / "heatmap_significant.png", dpi=600, bbox_inches="tight")
g.savefig(HM / "heatmap_significant.pdf", bbox_inches="tight")
g.savefig(HM / "heatmap_significant.svg", bbox_inches="tight")
plt.close(g.figure)

# C) Targeted heatmap of the spec-required genes
target_sym_set = TARGET_GENE_SYMBOLS
target_mat_idx = []
for sym in target_sym_set:
    ens_list = sym_to_ens.get(sym, [])
    for e in ens_list:
        if e in expr.index:
            target_mat_idx.append(e)
mat = expr.loc[target_mat_idx][col_order]
# replace index labels with display labels
def disp(e):
    sym = gene_map.get(e, e)
    for cat, lab in TARGETS:
        if lab.split(" (")[0] == sym:
            return lab
    return sym
new_idx = [disp(e) for e in mat.index]
mat.index = new_idx
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)
g = sns.clustermap(z, col_cluster=False, row_cluster=False, cmap="RdBu_r", center=0,
                   figsize=(7.5, 4.5), xticklabels=True, yticklabels=True,
                   col_colors=[col_colors_group, col_colors_patient, col_colors_batch])
g.ax_heatmap.set_title("Targeted genes — PDO vs PDC (z-score)", pad=14)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("target gene")
g.figure.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
g.savefig(HM / "heatmap_targeted.png", dpi=600, bbox_inches="tight")
g.savefig(HM / "heatmap_targeted.pdf", bbox_inches="tight")
g.savefig(HM / "heatmap_targeted.svg", bbox_inches="tight")
plt.close(g.figure)
# data csv for the targeted heatmap
mat.to_csv(HM / "heatmap_targeted.csv")

# ----- SINGLE-GENE BAR+SCATTER+PAIR FIGURES -----
def single_gene_figure(display_label, category, symbol):
    # collect per-sample normalized expression (CPM)
    ens_list = sym_to_ens.get(symbol, [])
    if not ens_list:
        # not detected
        fig, ax = plt.subplots(figsize=(5.5, 3.4))
        ax.text(0.5, 0.5, f"{display_label} not detected in filtered matrix",
                ha="center", va="center", fontsize=12)
        ax.axis("off")
        fig.tight_layout()
        for ext in ("png","pdf","svg"):
            fig.savefig(GE / f"{category}_{symbol}_not_detected.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
        return None
    ens = ens_list[0]
    if ens not in filt.index:
        return None
    # CPM
    libsize = filt.sum(axis=0)
    cpm_row = (filt.loc[ens] / libsize) * 1e6
    # rename to final names
    vals = pd.Series({meta.set_index("sample").loc[s, "final_sample_name"]: cpm_row[s]
                      for s in meta["sample"]})
    # group label
    groups = pd.Series({meta.set_index("sample").loc[s, "final_sample_name"]:
                        meta.set_index("sample").loc[s, "corrected_biology"]
                        for s in meta["sample"]})
    pairs = pd.Series({meta.set_index("sample").loc[s, "final_sample_name"]:
                       int(meta.set_index("sample").loc[s, "patient"])
                       for s in meta["sample"]})

    pdc_vals = vals[groups == "PDC"]
    pdo_vals = vals[groups == "PDO"]
    # Paired t-test (since PDC and PDO are patient-paired)
    common_patients = sorted(set(pairs.loc[pdc_vals.index]).intersection(pairs.loc[pdo_vals.index]))
    pdc_arr = np.array([pdc_vals.loc[[i for i in pdc_vals.index if pairs[i] == p]] for p in common_patients]).flatten()
    pdo_arr = np.array([pdo_vals.loc[[i for i in pdo_vals.index if pairs[i] == p]] for p in common_patients]).flatten()
    # because there's one value per patient per group
    pdc_per_patient = [float(pdc_vals[[i for i in pdc_vals.index if pairs[i]==p][0]]) for p in common_patients]
    pdo_per_patient = [float(pdo_vals[[i for i in pdo_vals.index if pairs[i]==p][0]]) for p in common_patients]
    if len(pdc_per_patient) >= 2:
        t_stat, p_paired = stats.ttest_rel(pdo_per_patient, pdc_per_patient)
        try:
            w_stat, p_wilcox = stats.wilcoxon(pdo_per_patient, pdc_per_patient)
        except Exception:
            w_stat, p_wilcox = np.nan, np.nan
        # use paired t-test as primary
        p_used = p_paired
        method = "paired t-test (PDO vs PDC, patient-paired)"
    else:
        p_used = np.nan; method = "insufficient pairs"
        w_stat = np.nan
    # also report FDR via BH over all targeted genes
    p_used = float(p_used) if not np.isnan(p_used) else np.nan
    # CI for diff in means
    diff = np.array(pdo_per_patient) - np.array(pdc_per_patient)
    if len(diff) >= 2:
        ci_low, ci_high = stats.t.interval(0.95, len(diff)-1, loc=diff.mean(), scale=stats.sem(diff))
    else:
        ci_low = ci_high = np.nan

    # Figure
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    bars = pd.DataFrame({
        "group": ["PDC", "PDO"],
        "mean":  [pdc_vals.mean(), pdo_vals.mean()],
        "sem":   [pdc_vals.std(ddof=1)/np.sqrt(len(pdc_vals)) if len(pdc_vals)>1 else 0,
                  pdo_vals.std(ddof=1)/np.sqrt(len(pdo_vals)) if len(pdo_vals)>1 else 0],
    })
    # bar with 95% CI
    n_pdc = len(pdc_vals); n_pdo = len(pdo_vals)
    def ci95(x):
        if len(x) < 2: return (np.nan, np.nan)
        return stats.t.interval(0.95, len(x)-1, loc=np.mean(x), scale=stats.sem(x))
    pdc_ci = ci95(pdc_vals.values); pdo_ci = ci95(pdo_vals.values)
    pdc_err = (pdc_vals.mean()-pdc_ci[0], pdc_ci[1]-pdc_vals.mean())
    pdo_err = (pdo_vals.mean()-pdo_ci[0], pdo_ci[1]-pdo_vals.mean())

    ax.bar("PDC", pdc_vals.mean(), color=PALETTE["PDC"], alpha=0.85,
           yerr=[[pdc_err[0]],[pdc_err[1]]], error_kw=dict(ecolor="black", lw=1.2, capsize=4))
    ax.bar("PDO", pdo_vals.mean(), color=PALETTE["PDO"], alpha=0.85,
           yerr=[[pdo_err[0]],[pdo_err[1]]], error_kw=dict(ecolor="black", lw=1.2, capsize=4))
    # scatter per-sample
    rng = np.random.default_rng(42)
    for grp, color in [("PDC", PALETTE["PDC"]), ("PDO", PALETTE["PDO"])]:
        idx = [c for c in vals.index if groups[c] == grp]
        x = np.full(len(idx), 0 if grp=="PDC" else 1)
        x = x + rng.uniform(-0.08, 0.08, len(idx))
        ax.scatter(x, vals[idx].values, s=40, color=color, edgecolor="black", zorder=3)
    # paired lines
    for p in common_patients:
        p_x0 = [i for i in pdc_vals.index if pairs[i]==p][0]
        p_x1 = [i for i in pdo_vals.index if pairs[i]==p][0]
        ax.plot([0, 1], [pdc_vals[p_x0], pdo_vals[p_x1]], "-", color="grey", lw=1.0, alpha=0.7, zorder=1)
    # significance bracket
    y_top = max(pdc_vals.max(), pdo_vals.max())
    y_range = y_top - min(pdc_vals.min(), pdo_vals.min())
    bracket_y = y_top + y_range*0.15
    ax.plot([0, 0, 1, 1], [bracket_y, bracket_y + y_range*0.05, bracket_y + y_range*0.05, bracket_y], color="black", lw=1.0)
    stars = ""
    if not np.isnan(p_used):
        if p_used < 0.0001: stars = "****"
        elif p_used < 0.001: stars = "***"
        elif p_used < 0.01:  stars = "**"
        elif p_used < 0.05:  stars = "*"
        else: stars = "ns"
    ax.text(0.5, bracket_y + y_range*0.07, stars, ha="center", va="bottom", fontsize=12)
    if not np.isnan(p_used):
        ax.text(0.5, bracket_y - y_range*0.05, f"p={p_used:.3g}", ha="center", va="top", fontsize=9, color="black")
    ax.set_xticks([0,1]); ax.set_xticklabels(["PDC", "PDO"])
    ax.set_ylabel("CPM")
    ax.set_title(f"{display_label}  ({category})\nPDO vs PDC")
    # caption text below plot
    caption = f"n_PDC={n_pdc}, n_PDO={n_pdo}, paired n={len(common_patients)}\n{method}; p={p_used if not np.isnan(p_used) else float('nan'):.3g}; diff 95% CI [{ci_low:.2g},{ci_high:.2g}] (CPM)"
    fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=8)
    fig.tight_layout()
    safe = f"{category}_{symbol}".replace(" ", "_").replace("(", "").replace(")", "").replace("/","_")
    for ext in ("png","pdf","svg"):
        fig.savefig(GE / f"{safe}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    # CSV
    pd.DataFrame({
        "sample": vals.index,
        "group": [groups[c] for c in vals.index],
        "patient": [int(pairs[c]) for c in vals.index],
        "CPM": vals.values,
    }).to_csv(GE / f"{safe}.csv", index=False)
    return {"p_paired": p_used, "method": method, "diff_CI": (ci_low, ci_high)}

gene_results = []
for (cat, lab), lookup_sym in zip(TARGETS, TARGET_LOOKUP_SYMBOLS):
    info = single_gene_figure(lab, cat, lookup_sym)
    gene_results.append({
        "category": cat, "label": lab, "gene_symbol": lookup_sym,
        "p_value": info["p_paired"] if info else np.nan,
        "method": info["method"] if info else "n/a",
        "diff_CI_low": info["diff_CI"][0] if info else np.nan,
        "diff_CI_high": info["diff_CI"][1] if info else np.nan,
    })
gd = pd.DataFrame(gene_results)
# BH-FDR across the gene-level p-values
from statsmodels.stats.multitest import multipletests
mask = gd["p_value"].notna()
gd["FDR"] = np.nan
gd.loc[mask, "FDR"] = multipletests(gd.loc[mask, "p_value"], method="fdr_bh")[1]
gd.to_csv(DE / "target_genes_single_sample_stats.csv", index=False)
gd.to_csv(DE / "target_genes_single_sample_stats.tsv", sep="\t", index=False)
print("\nSingle-sample paired stats:\n", gd.to_string(index=False))

# Save
print("STEP 3 complete")