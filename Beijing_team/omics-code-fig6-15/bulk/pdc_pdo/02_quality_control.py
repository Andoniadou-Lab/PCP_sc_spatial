"""Step 2 — Quality control + global structure figures (rewritten).

Inputs:
  - counts (raw, 62754 genes × 8 samples) from upstream featureCounts
  - sample metadata with corrected biology labels (PDC/PDO)

Outputs (under 01_quality_control/):
  qc_metrics.tsv
  sample_correlation_matrix.tsv / .csv
  pca_coordinates.tsv / .csv
  filtered_counts.tsv
  vst.tsv                          (log2(CPM+1) surrogate, used for QC plots & heatmaps)
  library_and_detected_genes.{png,pdf,svg,csv}
  correlation_heatmap.{png,pdf,svg}
  pca_group.{png,pdf,svg}
  pca_patient.{png,pdf,svg}
  pca_batch.{png,pdf,svg}
  hierarchical_clustering.{png,pdf,svg}
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA

TURN = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
QC = TURN / "01_quality_control"
QC.mkdir(parents=True, exist_ok=True)

UPSTREAM_META = "/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/count_matrix.metadata.csv"
UPSTREAM_RAW = "/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/count_matrix.raw.tsv"
GENE_INFO = "/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/gene_info.tsv"

RENAME = json.load(open(TURN / "00_metadata_audit" / "rename_dict.json"))

meta = pd.read_csv(UPSTREAM_META)
meta["final_sample_name"] = meta["sample"].map(RENAME)
counts = pd.read_csv(UPSTREAM_RAW, sep="\t", index_col=0)
gene_info = pd.read_csv(GENE_INFO, sep="\t")

# Trim version suffix
counts.index = [g.split(".")[0] for g in counts.index]
gene_info["ensembl"] = gene_info["Geneid"].str.split(".").str[0]
gene_info = gene_info.drop_duplicates("ensembl").set_index("ensembl")

# Filter: CPM >= 1 in >= 4 samples
libsize = counts.sum(axis=0)
cpm = counts.divide(libsize, axis=1) * 1e6
keep = (cpm >= 1).sum(axis=1) >= 4
filt = counts.loc[keep].copy()
filt = filt[meta["sample"].tolist()]
print(f"Genes: {counts.shape[0]} -> filtered {filt.shape[0]} (CPM>=1 in >=4 samples)")

# Per-sample QC
qc_rows = []
for s in meta["sample"]:
    sub = filt[s]
    detected = int((sub > 0).sum())
    qc_rows.append({
        "sample": s,
        "final_sample_name": meta.set_index("sample").loc[s, "final_sample_name"],
        "biology_corrected": "PDC" if meta.set_index("sample").loc[s, "biology"] == "spheres" else "PDO",
        "patient": int(meta.set_index("sample").loc[s, "patient"]),
        "batch": meta.set_index("sample").loc[s, "batch"],
        "library_size": int(libsize[s]),
        "detected_genes": detected,
        "median_counts_per_detected_gene": float(sub[sub > 0].median()) if detected else 0.0,
    })
qc = pd.DataFrame(qc_rows)
qc.to_csv(QC / "qc_metrics.tsv", sep="\t", index=False)

# log2(CPM+1) for QC/clustering/PCA/heatmap
libsize_vec = filt.sum(axis=0)
norm = filt / libsize_vec.values * 1e6
# Drop genes with NaN (e.g. unmapped IDs from upstream featureCounts)
norm = norm.dropna(how="any").fillna(0.0)
new_cols = [meta.set_index("sample").loc[c, "final_sample_name"] for c in norm.columns]
expr = pd.DataFrame(
    np.log2(norm.values + 1.0),
    index=norm.index,
    columns=new_cols,
).astype(float)
assert not expr.isnull().any().any(), "NaN in expr"
expr.to_csv(QC / "vst.tsv", sep="\t")
filt.loc[norm.index].to_csv(QC / "filtered_counts.tsv", sep="\t")

# === FIGURES ===
sns.set_theme(style="white", context="paper", font_scale=1.1)
PALETTE = {"PDC": "#3B82C4", "PDO": "#E07A5F"}

def save_fig(fig, base, subdir=QC):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(subdir / f"{base}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)

order = qc.sort_values("final_sample_name").reset_index(drop=True)

# 1. Library size + detected genes
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
sns.barplot(data=order, x="final_sample_name", y="library_size", hue="biology_corrected",
            palette=PALETTE, dodge=False, ax=axes[0])
axes[0].set_title("Library size per sample")
axes[0].set_ylabel("Total reads (counts)")
axes[0].set_xlabel("")
axes[0].tick_params(axis='x', rotation=45)
for i, r in order.iterrows():
    axes[0].text(i, r["library_size"] * 1.02, f'{r["library_size"]/1e6:.1f}M',
                  ha="center", fontsize=8)
axes[0].legend(title="")
sns.barplot(data=order, x="final_sample_name", y="detected_genes", hue="biology_corrected",
            palette=PALETTE, dodge=False, ax=axes[1])
axes[1].set_title("Detected genes per sample (counts > 0)")
axes[1].set_ylabel("# detected genes")
axes[1].set_xlabel("")
axes[1].tick_params(axis='x', rotation=45)
axes[1].legend(title="")
fig.tight_layout()
save_fig(fig, "library_and_detected_genes")
order[["final_sample_name", "biology_corrected", "library_size", "detected_genes"]].to_csv(
    QC / "library_and_detected_genes.csv", index=False)

# 2. Sample correlation heatmap
corr = expr.corr(method="pearson")
corr.to_csv(QC / "sample_correlation_matrix.tsv", sep="\t")
corr.to_csv(QC / "sample_correlation_matrix.csv")
fig, ax = plt.subplots(figsize=(6, 5))
sns.heatmap(corr, annot=True, fmt=".3f", cmap="viridis", vmin=corr.values.min(), vmax=1.0,
            square=True, cbar_kws={"label": "Pearson r"}, ax=ax, linewidths=0.5)
ax.set_title("Sample correlation (Pearson, log2(CPM+1))")
fig.tight_layout()
save_fig(fig, "correlation_heatmap")

# 3. PCA
pca = PCA(n_components=4)
pcs = pca.fit_transform(expr.T.values)
pca_df = pd.DataFrame(pcs[:, :4], index=expr.columns,
                      columns=["PC1", "PC2", "PC3", "PC4"])
pca_df["biology"] = [qc.set_index("final_sample_name").loc[s, "biology_corrected"] for s in pca_df.index]
pca_df["patient"] = [qc.set_index("final_sample_name").loc[s, "patient"] for s in pca_df.index]
pca_df["batch"] = [qc.set_index("final_sample_name").loc[s, "batch"] for s in pca_df.index]
pca_df.to_csv(QC / "pca_coordinates.tsv", sep="\t")
pca_df.to_csv(QC / "pca_coordinates.csv")

fig, ax = plt.subplots(figsize=(6.5, 5))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="biology", palette=PALETTE,
                s=180, ax=ax, edgecolor="black")
for s, r in pca_df.iterrows():
    ax.text(r["PC1"]+0.5, r["PC2"]+0.5, s, fontsize=9)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA — PDO vs PDC (log2(CPM+1))")
ax.legend(title="")
fig.tight_layout()
save_fig(fig, "pca_group")

fig, ax = plt.subplots(figsize=(6.5, 5))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="patient", palette="Set2",
                s=180, ax=ax, edgecolor="black")
for s, r in pca_df.iterrows():
    ax.text(r["PC1"]+0.5, r["PC2"]+0.5, s, fontsize=9)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA — by patient")
ax.legend(title="patient")
fig.tight_layout()
save_fig(fig, "pca_patient")

fig, ax = plt.subplots(figsize=(6.5, 5))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="batch", palette="Set1",
                s=180, ax=ax, edgecolor="black")
for s, r in pca_df.iterrows():
    ax.text(r["PC1"]+0.5, r["PC2"]+0.5, s, fontsize=9)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA — by sequencing batch")
ax.legend(title="batch")
fig.tight_layout()
save_fig(fig, "pca_batch")

# 4. Hierarchical clustering
fig, ax = plt.subplots(figsize=(7, 4))
Z = linkage(expr.T.values, method="average", metric="euclidean")
dendrogram(Z, labels=list(expr.columns), ax=ax, leaf_font_size=10, above_threshold_color="grey")
ax.set_title("Hierarchical clustering — samples (Euclidean, average linkage)")
ax.set_ylabel("Distance")
fig.tight_layout()
save_fig(fig, "hierarchical_clustering")

# Save augmented metadata for downstream
meta.to_csv(TURN / "00_metadata_audit" / "sample_metadata_with_final_names.tsv", sep="\t", index=False)

print("STEP 2 (QC) complete")
print(qc.to_string(index=False))