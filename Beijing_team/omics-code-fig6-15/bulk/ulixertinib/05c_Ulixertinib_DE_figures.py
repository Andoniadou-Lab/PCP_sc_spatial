"""Step 5d (part 1) — Ulixertinib vs Ctrl: featureCounts -> count matrix ->
QC -> DESeq2 -> volcano -> heatmaps -> single-gene.

Runs AFTER 06_scripts/featureCounts_ulix.sh finishes all 8 BAM samples.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime, timezone
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

TURN = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
A2 = TURN / "03_Ulixertinib_vs_Ctrl"
FC_DIR = A2 / "featureCounts"

DE = A2 / "differential_expression"; DE.mkdir(exist_ok=True)
QC = A2 / "quality_control"; QC.mkdir(exist_ok=True)
HM = A2 / "heatmap"; HM.mkdir(exist_ok=True)
GS = A2 / "gene_expression"; GS.mkdir(exist_ok=True)
GE = A2 / "GSEA"; GE.mkdir(exist_ok=True)

TARGET_GENES = ["SOX2","SOX9","FLRT2","STOX2","MKI67","TACSTD2","CASP3","BCL2"]
TARGET_LABELS = ["SOX2","SOX9","FLRT2","STOX2","MKI67","TROP2 (TACSTD2)","CASP3","BCL2"]
TARGET_CAT = {
    "SOX2":"Stemness","SOX9":"Stemness","FLRT2":"Stemness","STOX2":"Stemness",
    "MKI67":"Proliferation","TACSTD2":"Proliferation","CASP3":"Apoptosis","BCL2":"Apoptosis",
}
def cat_for(sym): return TARGET_CAT.get(sym, "Other")

# ---- 1. Build count matrix ----
sheet = pd.read_csv(A2 / "sample_sheet.csv")
samples = sheet["sample"].tolist()
bam_map = dict(zip(sheet["sample"], sheet["bam"]))

counts = None
ref_genes = None
for s in samples:
    bam = bam_map[s]
    fc_path = FC_DIR / f"{bam.replace('.bam','')}.featureCounts.tsv"
    if not fc_path.exists():
        raise FileNotFoundError(fc_path)
    df = pd.read_csv(fc_path, sep="\t", comment="#")
    if counts is None:
        counts = pd.DataFrame(index=df["Geneid"].values)
        ref_genes = df[["Geneid","gene_name"]].rename(columns={"Geneid":"ensembl","gene_name":"gene_symbol"})
    counts[s] = df[df.columns[-1]].values

print(f"count matrix: {counts.shape[0]} x {counts.shape[1]}")
counts.index.name = "ensembl"
counts.to_csv(A2 / "count_matrix.raw.tsv", sep="\t")
ref_genes.to_csv(A2 / "count_matrix.metadata.tsv", sep="\t", index=False)
counts.T.to_csv(A2 / "count_matrix.csv")

# ---- 2. Filter ----
libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
keep = (cpm >= 1).sum(axis=1) >= 4
filtered = counts.loc[keep].copy()
print(f"after filter: {filtered.shape[0]} genes")
filtered.to_csv(A2 / "count_matrix.filtered.tsv", sep="\t")
log_cpm = np.log2(cpm.loc[keep] + 1)
log_cpm.T.to_csv(A2 / "vst.tsv", sep="\t")

# ---- 3. QC ----
sns.set_style("whitegrid")
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
colors = ["#4DBBD5" if g == "Ulixertinib" else "#E64B35" for g in sheet["group"]]
ax[0].bar(samples, libsize/1e6, color=colors); ax[0].set_ylabel("Library size (M read pairs)"); ax[0].set_title("Library size"); ax[0].tick_params(axis="x", rotation=45)
n_det = (filtered > 0).sum(axis=0)
ax[1].bar(samples, n_det, color=colors); ax[1].set_ylabel("Genes detected (filtered)"); ax[1].set_title("Genes detected"); ax[1].tick_params(axis="x", rotation=45)
fig.suptitle("Ulixertinib vs Ctrl — library size & detected genes")
fig.tight_layout()
for ext in ("png","pdf","svg"): fig.savefig(QC / f"library_and_detected_genes.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)
pd.DataFrame({"sample":samples,"group":sheet["group"],"libsize":libsize.values,"n_detected":n_det.values}).to_csv(QC / "library_and_detected_genes.csv", index=False)

corr = log_cpm.T.corr()
fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(corr, annot=True, fmt=".2f", cmap="viridis", vmin=0.8, vmax=1.0, ax=ax)
ax.set_title("Sample correlation — Ulixertinib vs Ctrl (log2(CPM+1))")
fig.tight_layout()
for ext in ("png","pdf","svg"): fig.savefig(QC / f"correlation_heatmap.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)
corr.to_csv(QC / "sample_correlation_matrix.csv")

pca = PCA(n_components=4)
pcs = pca.fit_transform(log_cpm.T.values)
pca_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(4)])
pca_df["sample"] = log_cpm.columns
pca_df["group"] = sheet["group"].values
pca_df["patient"] = sheet["patient_id"].values
pca_df.to_csv(QC / "pca_coordinates.csv", index=False)

for by, title in [("group","Group (Ulixertinib vs Ctrl)"),("patient","Patient")]:
    fig, ax = plt.subplots(figsize=(6.5, 5))
    if by == "group":
        for g, sub in pca_df.groupby("group"):
            ax.scatter(sub["PC1"], sub["PC2"], label=g, s=140, edgecolor="k",
                       c="#4DBBD5" if g=="Ulixertinib" else "#E64B35")
    else:
        for p, sub in pca_df.groupby("patient"):
            ax.scatter(sub["PC1"], sub["PC2"], label=f"Patient {p}", s=140, edgecolor="k")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title(f"PCA coloured by {title}")
    ax.legend()
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(QC / f"pca_{by}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
print("QC done")

# ---- 4. DESeq2 paired design ----
meta = sheet.set_index("sample")[["group","patient_id"]].rename(columns={"patient_id":"patient"})
meta["group"] = pd.Categorical(meta["group"], categories=["Ctrl","Ulixertinib"])
mat_for_de = filtered.T
dds = DeseqDataSet(counts=mat_for_de, metadata=meta, design_factors=["patient","group"], refit_cooks=False, n_cpus=1)
dds.deseq2()
stat = DeseqStats(dds, contrast=["group","Ulixertinib","Ctrl"], n_cpus=1)
stat.summary()
res = stat.results_df.copy().reset_index().rename(columns={"index":"ensembl"})
res = res.merge(ref_genes, on="ensembl", how="left")

def direction(r):
    if pd.isna(r["padj"]): return "ns"
    if r["padj"] < 0.05 and abs(r["log2FoldChange"]) >= 1:
        return "up_in_Ulixertinib" if r["log2FoldChange"] > 0 else "up_in_Ctrl"
    return "ns"
res["direction"] = res.apply(direction, axis=1)
res["significance"] = res["direction"].apply(lambda d: "FDR<0.05 & |log2FC|>=1" if d != "ns" else "ns")
res.to_csv(DE / "all_genes_Ulixertinib_vs_Ctrl.tsv", sep="\t", index=False)
sig = res[res["direction"] != "ns"].copy()
sig.to_csv(DE / "significant_genes_Ulixertinib_vs_Ctrl.tsv", sep="\t", index=False)
print(f"DE: up_in_Ulixertinib={(res.direction=='up_in_Ulixertinib').sum()}, up_in_Ctrl={(res.direction=='up_in_Ctrl').sum()}, ns={(res.direction=='ns').sum()}")

ranked = res.dropna(subset=["stat"]).sort_values("stat", ascending=False)
ranked[["gene_symbol","stat"]].rename(columns={"gene_symbol":"#Name","stat":"score"}).to_csv(GE / "prerank_input.tsv", sep="\t", index=False)

# ---- 5. Volcano ----
v_df = res.dropna(subset=["log2FoldChange","padj"]).copy()
v_df["nlog10_p"] = -np.log10(v_df["padj"].clip(lower=1e-300))
v_colors = []
for _, r in v_df.iterrows():
    if r["padj"] < 0.05 and r["log2FoldChange"] >= 1: v_colors.append("#D6604D")
    elif r["padj"] < 0.05 and r["log2FoldChange"] <= -1: v_colors.append("#4393C3")
    else: v_colors.append("#BBBBBB")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(v_df["log2FoldChange"], v_df["nlog10_p"], c=v_colors, s=14, alpha=0.7, edgecolor="none")
ax.axhline(-np.log10(0.05), ls="--", color="grey", lw=0.8)
ax.axvline(1, ls="--", color="grey", lw=0.8); ax.axvline(-1, ls="--", color="grey", lw=0.8)
for sym, lab in zip(TARGET_GENES, TARGET_LABELS):
    sub = v_df[v_df["gene_symbol"] == sym]
    if not sub.empty:
        for _, r in sub.iterrows():
            ax.text(r["log2FoldChange"], r["nlog10_p"], f" {lab}", fontsize=9, fontweight="bold")
ax.set_xlabel("log2 fold change (Ulixertinib vs Ctrl)")
ax.set_ylabel("-log10(FDR)")
up_n = (v_df["direction"]=="up_in_Ulixertinib").sum()
dn_n = (v_df["direction"]=="up_in_Ctrl").sum()
ax.set_title(f"Ulixertinib vs Ctrl — up_in_Ulixertinib={up_n}, up_in_Ctrl={dn_n}")
fig.tight_layout()
for ext in ("png","pdf","svg"): fig.savefig(DE / f"volcano_Ulixertinib_vs_Ctrl.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)
v_df.to_csv(DE / "volcano_Ulixertinib_vs_Ctrl.csv", index=False)

# ---- 6. Heatmaps ----
def heatmap(df_z, name, title):
    fig, ax = plt.subplots(figsize=(8, max(6, 0.18*df_z.shape[0])))
    sns.heatmap(df_z, cmap="RdBu_r", center=0, ax=ax, yticklabels=True if df_z.shape[0]<=60 else False)
    col_colors = ["#E64B35" if g=="Ctrl" else "#4DBBD5" for g in sheet["group"]]
    for i, c in enumerate(col_colors):
        ax.add_patch(plt.Rectangle((i, -0.5), 1, 0.3, color=c, transform=ax.transData, clip_on=False))
    ax.set_title(title)
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(HM / f"{name}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    df_z.to_csv(HM / f"{name}.csv")

z = log_cpm.sub(log_cpm.mean(axis=1), axis=0).div(log_cpm.std(axis=1).replace(0,1), axis=0)
top_var = z.var(axis=1).sort_values(ascending=False).head(1000).index
heatmap(z.loc[top_var].T[samples].T, "heatmap_full", "Top 1000 variable genes — Ulixertinib vs Ctrl")
sig_idx = sig["ensembl"].values
heatmap(z.loc[z.index.isin(sig_idx)].T[samples].T, "heatmap_significant", "Significant DE genes — Ulixertinib vs Ctrl")
target_idx = ref_genes[ref_genes["gene_symbol"].isin(TARGET_GENES)]["ensembl"]
heatmap(z.loc[z.index.isin(target_idx)].T[samples].T, "heatmap_targeted", "Targeted genes — Ulixertinib vs Ctrl")

# ---- 7. Single-gene plots ----
def single_gene(symbol, display):
    ens_row = ref_genes[ref_genes["gene_symbol"] == symbol]
    if ens_row.empty: return
    ens = ens_row["ensembl"].iloc[0]
    if ens not in counts.index: return
    vals_full = (counts.loc[ens] / counts.sum(axis=0) * 1e6).reindex(samples).values
    if (vals_full < 1).all():
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, f"{symbol}\n(low expression in all 8 samples — CPM<1)", ha="center", va="center", fontsize=11)
        ax.axis("off")
        safe = display.replace(" ","_").replace("(","").replace(")","")
        for ext in ("png","pdf","svg"): fig.savefig(GS / f"{safe}_low_expression.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
        return
    grp = sheet["group"].values
    pat = sheet["patient_id"].values
    pat_levels = sorted(pat[grp=="Ctrl"].tolist())
    cp, up = [], []
    for p in pat_levels:
        cp.append(float(vals_full[(grp=="Ctrl") & (pat==p)][0]))
        up.append(float(vals_full[(grp=="Ulixertinib") & (pat==p)][0]))
    diffs = np.array(up) - np.array(cp)
    t, p = stats.ttest_rel(up, cp)
    n = len(diffs); sd = diffs.std(ddof=1); sem = sd/np.sqrt(n)
    ci_l, ci_h = stats.t.interval(0.95, n-1, loc=diffs.mean(), scale=sem)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    means = [np.mean(cp), np.mean(up)]
    sems = [np.std(cp, ddof=1)/np.sqrt(len(cp)), np.std(up, ddof=1)/np.sqrt(len(up))]
    ax.bar([0,1], means, yerr=sems, color=["#E64B35","#4DBBD5"], capsize=8, edgecolor="k")
    for i in range(len(pat_levels)):
        ax.scatter([0],[cp[i]],color="white",edgecolor="k",s=60,zorder=3)
        ax.scatter([1],[up[i]],color="white",edgecolor="k",s=60,zorder=3)
        ax.plot([0,1],[cp[i], up[i]], color="grey", lw=0.7, zorder=2)
    sig_lab = "ns"
    if p < 0.0001: sig_lab = "****"
    elif p < 0.001: sig_lab = "***"
    elif p < 0.01: sig_lab = "**"
    elif p < 0.05: sig_lab = "*"
    ymax = max(max(cp), max(up)) * 1.30
    ax.plot([0,1],[ymax*0.95, ymax*0.95], color="k", lw=1)
    ax.text(0.5, ymax*1.0, f"{sig_lab}\npaired t p={p:.3g}\nΔ={diffs.mean():.2f} [{ci_l:.2f},{ci_h:.2f}]",
            ha="center", va="bottom", fontsize=9)
    ax.set_xticks([0,1]); ax.set_xticklabels(["Ctrl (n=4)","Ulixertinib (n=4)"])
    ax.set_ylabel(f"{symbol} CPM")
    ax.set_title(f"{display} — Ulixertinib vs Ctrl (paired by patient)")
    ax.set_ylim(0, ymax*1.25)
    fig.tight_layout()
    safe = display.replace(" ","_").replace("(","").replace(")","")
    for ext in ("png","pdf","svg"): fig.savefig(GS / f"{safe}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame({
        "patient":pat_levels, "Ctrl":cp, "Ulixertinib":up, "diff":diffs,
        "p_paired_t":p, "mean_diff":diffs.mean(), "CI_low":ci_l, "CI_high":ci_h,
    }).to_csv(GS / f"{safe}.csv", index=False)

target_rows = []
for sym, lab in zip(TARGET_GENES, TARGET_LABELS):
    r = res[res["gene_symbol"] == sym]
    if r.empty:
        target_rows.append({"category":cat_for(sym),"display_label":lab,"gene_symbol":sym,
                            "ensembl":None,"log2FC":None,"SE":None,"Wald_statistic":None,
                            "P_value":None,"FDR":None,"baseMean":None,
                            "detection_status":"filtered_due_to_low_expression","low_expression_filtered":"yes",
                            "duplicate_gene_symbol":"no"})
        continue
    rr = r.iloc[0]
    det = "yes" if rr["baseMean"] > 0 else "filtered_due_to_low_expression"
    target_rows.append({"category":cat_for(sym),"display_label":lab,"gene_symbol":sym,
                        "ensembl":rr["ensembl"],"log2FC":rr["log2FoldChange"],"SE":rr["lfcSE"],
                        "Wald_statistic":rr["stat"],"P_value":rr["pvalue"],"FDR":rr["padj"],
                        "baseMean":rr["baseMean"],"detection_status":det,
                        "low_expression_filtered":"yes" if det=="filtered_due_to_low_expression" else "no",
                        "duplicate_gene_symbol":"no"})
pd.DataFrame(target_rows).to_csv(DE / "target_genes_Ulixertinib_vs_Ctrl.tsv", sep="\t", index=False)

single_rows = []
for sym in TARGET_GENES:
    ens_row = ref_genes[ref_genes["gene_symbol"] == sym]
    if ens_row.empty: continue
    ens = ens_row["ensembl"].iloc[0]
    if ens not in counts.index: continue
    vals = (counts.loc[ens] / counts.sum(axis=0) * 1e6).reindex(samples).values
    grp = sheet["group"].values
    pat = sheet["patient_id"].values
    if (vals < 1).all():
        single_rows.append({"label":sym,"gene_symbol":sym,"category":cat_for(sym),
                            "p_value":np.nan,"method":"paired t-test (Ulixertinib vs Ctrl, patient-paired)",
                            "diff_CI_low":np.nan,"diff_CI_high":np.nan})
        continue
    pat_lv = sorted(pat[grp=="Ctrl"].tolist())
    cp, up = [], []
    for p in pat_lv:
        cp.append(float(vals[(grp=="Ctrl")&(pat==p)][0]))
        up.append(float(vals[(grp=="Ulixertinib")&(pat==p)][0]))
    diffs = np.array(up) - np.array(cp)
    t, pv = stats.ttest_rel(up, cp)
    n = len(diffs); sd = diffs.std(ddof=1); sem = sd/np.sqrt(n)
    ci_l, ci_h = stats.t.interval(0.95, n-1, loc=diffs.mean(), scale=sem)
    single_rows.append({"label":sym,"gene_symbol":sym,"category":cat_for(sym),
                        "p_value":pv,"method":"paired t-test (Ulixertinib vs Ctrl, patient-paired)",
                        "diff_CI_low":ci_l,"diff_CI_high":ci_h})

sg_df = pd.DataFrame(single_rows)
if not sg_df.empty and sg_df["p_value"].notna().any():
    pvals = sg_df["p_value"].fillna(1.0).values
    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    sg_df["FDR"] = fdr
sg_df.to_csv(DE / "target_genes_single_sample_stats.csv", index=False)

for sym, lab in zip(TARGET_GENES, TARGET_LABELS):
    single_gene(sym, lab)

print("DONE step 5d part 1")