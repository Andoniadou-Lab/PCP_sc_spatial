"""Step 4 — Analysis 1: pathway scores, GSEA, GO/KEGG for PDO vs PDC."""
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
from statsmodels.stats.multitest import multipletests

TURN = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
QC = TURN / "01_quality_control"
A1 = TURN / "02_PDC_vs_PDO"
DE = A1 / "differential_expression"
PS = A1 / "pathway_scores"
GSEA = A1 / "GSEA"
GK = A1 / "GO_KEGG"
for p in [PS, GSEA, GK]: p.mkdir(parents=True, exist_ok=True)

RENAME = json.load(open(TURN / "00_metadata_audit" / "rename_dict.json"))
meta = pd.read_csv(TURN / "00_metadata_audit" / "sample_metadata_with_final_names.tsv", sep="\t")
meta["corrected_biology"] = meta["biology"].map({"single_cells": "PDO", "spheres": "PDC"})
expr = pd.read_csv(QC / "vst.tsv", sep="\t", index_col=0)
res_df = pd.read_csv(DE / "all_genes_PDO_vs_PDC.tsv", sep="\t")
gene_info = pd.read_csv("/home/susu/Rprojects/PCP202607/outputs/turn_4e5b864f455947418f606e79d0548cf7/featureCounts/gene_info.tsv", sep="\t")
gene_info["ensembl"] = gene_info["Geneid"].str.split(".").str[0]
gene_info = gene_info.drop_duplicates("ensembl").set_index("ensembl")
gene_map = gene_info["gene_name"].to_dict()
PALETTE = {"PDC": "#3B82C4", "PDO": "#E07A5F"}

# Curated gene sets (each line cites a public source)
GENE_SETS = {
    "Stemness (Wong 2008 + Bhatt 2012)": {
        "source": "Wong DJ 2008 Genes Dev; Bhatt K 2012 (ESC core + iPSC factors)",
        "genes": ["POU5F1","SOX2","NANOG","KLF4","MYC","LIN28A","LIN28B","NR0B1","ZFP42","TERT",
                  "UTF1","DPPA2","DPPA4","DPPA5","SALL4","ESRRB","ZIC3","TDGF1","FGF2","LEFTY1",
                  "LEFTY2","NODAL","FGF4","FGF5","GATA6","GATA4","TBX3","OTX2","PRDM14"],
    },
    "Proliferation (Hallmark E2F + G2M)": {
        "source": "Liberzon 2015 Hallmark E2F_TARGETS + G2M_CHECKPOINT (MSigDB v7)",
        "genes": ["MKI67","PCNA","CCNB1","CCNB2","CCND1","CCNE1","CCNE2","CCNA2","CDK1","CDK2","CDK4","CDK6",
                  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","AURKA","AURKB","BIRC5","TOP2A","CENPA","CENPB",
                  "BUB1","BUB1B","MAD2L1","TTK","KIF2C","KIF23","NDC80","SPC24","SPC25","UBE2C","UBE2S",
                  "CDC20","CDCA8","CDCA3","CDCA5","CDCA7","CDCA7L","ORC1","ORC6","GMNN","MCM10",
                  "POLD1","POLD2","POLE","POLE2","MRE11","RAD51","TYMS","RRM1","RRM2","FEN1","LIG1",
                  "TUBB","TUBA1B","STMN1","STMN2","KPNA2","XPO1","H2AX","H2AZ1","HMGB1","HMGB2"],
    },
    "Apoptosis (Hallmark Apoptosis core)": {
        "source": "Liberzon 2015 Hallmark APOPTOSIS (MSigDB v7)",
        "genes": ["CASP3","CASP6","CASP7","CASP8","CASP9","CASP10","BAX","BAK1","BCL2","BCL2L1","BCL2L11",
                  "BID","BAD","BIRC2","BIRC3","BIRC5","XIAP","CFLAR","FADD","TRADD","TNFRSF10A","TNFRSF10B",
                  "TNFRSF1A","FAS","DIABLO","CYCS","APAF1","AIFM1","ENDOG","PARP1","PARP2","CASP4",
                  "CASP5","CASP1","GSDMD","PYCARD","NLRP3","AIM2"],
    },
    "MAPK/ERK cascade (KEGG hsa04010)": {
        "source": "KEGG hsa04010 MAPK signaling; Reactome RAF/MAPK cascade",
        "genes": ["MAPK1","MAPK3","MAPK8","MAPK9","MAPK10","MAPK11","MAPK14","MAP2K1","MAP2K2","MAP2K7",
                  "MAP2K4","MAP3K1","MAP3K7","MAP3K11","RAF1","ARAF","BRAF","KRAS","NRAS","HRAS","SOS1","SOS2",
                  "GRB2","SHC1","PTPN11","RASGRP1","RASGRP2","RASGRP3","RAP1A","RAP1B","MKNK1","MKNK2",
                  "RPS6KA1","RPS6KA2","RPS6KA3","RPS6KA4","RPS6KA5","DUSP1","DUSP4","DUSP5","DUSP6","ELK1",
                  "ETS1","JUN","FOS","MYC","CCND1","EGFR","ERBB2","ERBB3","ERBB4","FGFR1","FGFR2","FGFR3",
                  "PDGFRA","PDGFRB","VEGFA","VEGFB","VEGFC","FGF2","FGF7","HGF","MET"],
    },
    "Cell cycle (GO 0007049)": {
        "source": "GO BP 0007049 cell cycle; Reactome Cell Cycle",
        "genes": ["CDK1","CDK2","CDK4","CDK6","CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2",
                  "CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","E2F1","E2F2","E2F3","E2F4","E2F5",
                  "RB1","RBL1","RBL2","TP53","MDM2","ATM","ATR","CHEK1","CHEK2","WEE1","CDC25A","CDC25B","CDC25C"],
    },
}

sym_to_ens = {}
for ens, sym in gene_map.items():
    sym_to_ens.setdefault(sym, []).append(ens)

# Per-sample mean-of-z score for each set
scores = {}
for set_name, info in GENE_SETS.items():
    ens_in = []
    for g in info["genes"]:
        for e in sym_to_ens.get(g, []):
            if e in expr.index: ens_in.append(e)
    if not ens_in:
        print("WARN: empty set:", set_name); continue
    scores[set_name] = expr.loc[ens_in].mean(axis=0)
    print(f"{set_name}: {len(ens_in)}/{len(info['genes'])} mapped")
scores_df = pd.DataFrame(scores)
scores_df.index.name = "final_sample_name"
scores_df.to_csv(PS / "pathway_scores.tsv", sep="\t")
scores_df.to_csv(PS / "pathway_scores.csv")

# paired test per pathway + bar+scatter+pair figure
score_results = []
for set_name in scores_df.columns:
    vals = scores_df[set_name]
    pdc_names = [c for c in vals.index if c.startswith("PDC")]
    pdo_names = [c for c in vals.index if c.startswith("PDO")]
    pairs_pdc = {int(meta.set_index("final_sample_name").loc[n, "patient"]): vals[n] for n in pdc_names}
    pairs_pdo = {int(meta.set_index("final_sample_name").loc[n, "patient"]): vals[n] for n in pdo_names}
    common = sorted(set(pairs_pdc).intersection(pairs_pdo))
    pdc_arr = np.array([pairs_pdc[p] for p in common])
    pdo_arr = np.array([pairs_pdo[p] for p in common])
    if len(pdc_arr) >= 2:
        t_stat, p_paired = stats.ttest_rel(pdo_arr, pdc_arr)
        diff = pdo_arr - pdc_arr
        ci_low, ci_high = stats.t.interval(0.95, len(diff)-1, loc=diff.mean(), scale=stats.sem(diff))
    else:
        p_paired = np.nan; ci_low = np.nan; ci_high = np.nan; diff = pdo_arr - pdc_arr
    score_results.append({"pathway": set_name, "n_PDC": len(pdc_arr), "n_PDO": len(pdo_arr),
                          "PDC_mean": float(pdc_arr.mean()), "PDO_mean": float(pdo_arr.mean()),
                          "PDO_minus_PDC_diff": float(diff.mean()),
                          "diff_CI_low": ci_low, "diff_CI_high": ci_high,
                          "paired_t_pvalue": p_paired, "genes_in_set": int(scores_df.shape[0])})
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    def ci(x):
        if len(x) < 2: return (np.nan, np.nan)
        return stats.t.interval(0.95, len(x)-1, loc=np.mean(x), scale=stats.sem(x))
    pdc_ci = ci(pdc_arr); pdo_ci = ci(pdo_arr)
    pdc_err = (pdc_arr.mean()-pdc_ci[0], pdc_ci[1]-pdc_arr.mean())
    pdo_err = (pdo_arr.mean()-pdo_ci[0], pdo_ci[1]-pdo_arr.mean())
    ax.bar("PDC", pdc_arr.mean(), color=PALETTE["PDC"], alpha=0.85, yerr=[[pdc_err[0]],[pdc_err[1]]], error_kw=dict(ecolor="black", lw=1.2, capsize=4))
    ax.bar("PDO", pdo_arr.mean(), color=PALETTE["PDO"], alpha=0.85, yerr=[[pdo_err[0]],[pdo_err[1]]], error_kw=dict(ecolor="black", lw=1.2, capsize=4))
    rng = np.random.default_rng(42)
    for grp, color, arr in [("PDC", PALETTE["PDC"], pdc_arr), ("PDO", PALETTE["PDO"], pdo_arr)]:
        x = np.full(len(arr), 0 if grp=="PDC" else 1) + rng.uniform(-0.08, 0.08, len(arr))
        ax.scatter(x, arr, s=45, color=color, edgecolor="black", zorder=3)
    for p in common:
        ax.plot([0, 1], [pairs_pdc[p], pairs_pdo[p]], "-", color="grey", lw=1.0, alpha=0.7, zorder=1)
    y_top = max(pdc_arr.max(), pdo_arr.max()); y_bot = min(pdc_arr.min(), pdo_arr.min())
    y_range = y_top - y_bot
    bracket_y = y_top + y_range*0.15
    ax.plot([0,0,1,1], [bracket_y, bracket_y+y_range*0.05, bracket_y+y_range*0.05, bracket_y], color="black", lw=1.0)
    stars = "ns"
    if not np.isnan(p_paired):
        if p_paired < 0.0001: stars = "****"
        elif p_paired < 0.001: stars = "***"
        elif p_paired < 0.01:  stars = "**"
        elif p_paired < 0.05:  stars = "*"
    ax.text(0.5, bracket_y+y_range*0.07, stars, ha="center", fontsize=12)
    if not np.isnan(p_paired):
        ax.text(0.5, bracket_y-y_range*0.05, f"p={p_paired:.3g}", ha="center", fontsize=9)
    ax.set_xticks([0,1]); ax.set_xticklabels(["PDC","PDO"])
    ax.set_ylabel("pathway score (mean log2(CPM+1))")
    ax.set_title(f"{set_name}\nPDO vs PDC (patient-paired)")
    fig.text(0.5, -0.02, f"n_PDC={len(pdc_arr)}, n_PDO={len(pdo_arr)}, paired n={len(common)}\npaired t-test; p={p_paired:.3g}; diff 95% CI [{ci_low:.2g},{ci_high:.2g}]",
             ha="center", fontsize=8)
    fig.tight_layout()
    safe = set_name.replace(" ","_").replace("(","").replace(")","").replace("+","_").replace("/","_")
    for ext in ("png","pdf","svg"): fig.savefig(PS / f"score_{safe}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame({"sample": list(pdc_names)+list(pdo_names),
                  "group": ["PDC"]*len(pdc_names)+["PDO"]*len(pdo_names),
                  "patient": [int(meta.set_index("final_sample_name").loc[s,"patient"]) for s in list(pdc_names)+list(pdo_names)],
                  "score": list(pdc_arr)+list(pdo_arr)}).to_csv(PS / f"score_{safe}.csv", index=False)

ps_df = pd.DataFrame(score_results)
mask = ps_df["paired_t_pvalue"].notna()
ps_df["FDR"] = np.nan
ps_df.loc[mask, "FDR"] = multipletests(ps_df.loc[mask, "paired_t_pvalue"], method="fdr_bh")[1]
ps_df.to_csv(PS / "pathway_scores_stats.tsv", sep="\t", index=False)
ps_df.to_csv(PS / "pathway_scores_stats.csv", index=False)
print("\nPathway score stats:\n", ps_df.to_string(index=False))

# ====== GSEA via gseapy prerank ======
import gseapy as gp
rnk = res_df.dropna(subset=["Wald_statistic"]).set_index("gene_symbol")["Wald_statistic"].sort_values(ascending=False)
rnk = rnk[~rnk.index.duplicated(keep="first")]

# Curated sets for prerank
gmt = {k: [g for g in v["genes"] if g in sym_to_ens and any(e in expr.index for e in sym_to_ens[g])] for k, v in GENE_SETS.items()}

# Try fetching MSigDB Hallmark via gseapy
HALLMARK = None
try:
    HALLMARK = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
except Exception as e:
    print("Hallmark fetch failed:", e)

def run_prerank(rnk, gene_sets_dict, tag):
    rows = []
    for k, genes in gene_sets_dict.items():
        if not genes: continue
        try:
            res = gp.prerank(rnk=rnk, gene_sets={k: genes}, outdir=None,
                             min_size=5, max_size=2000, permutation_num=1000,
                             seed=42, verbose=False, no_plot=True)
            r = res.res2d
            if r is None or r.empty: continue
            r["library"] = k
            rows.append(r)
        except Exception as e:
            print(f"prerank fail {k}: {e}")
    if not rows: return pd.DataFrame()
    full = pd.concat(rows, ignore_index=True)
    full.to_csv(GSEA / f"prerank_{tag}.tsv", sep="\t", index=False)
    return full

print("Running GSEA (curated)…")
gsea_cur = run_prerank(rnk, gmt, "curated")
print("Running GSEA (Hallmark)…")
gsea_hall = run_prerank(rnk, HALLMARK, "hallmark") if HALLMARK else pd.DataFrame()

all_gsea = pd.concat([gsea_cur, gsea_hall], ignore_index=True) if (not gsea_cur.empty or not gsea_hall.empty) else pd.DataFrame()
if not all_gsea.empty:
    all_gsea.to_csv(GSEA / "gsea_all_pathways.tsv", sep="\t", index=False)

# Single-pathway GSEA plots for focused sets
FOCUS = list(GENE_SETS.keys())
def gsea_plot(rnk_series, gene_set, set_name, out_base, nes, pval, fdr, leading_edge):
    genes_in_set = [g for g in gene_set if g in rnk_series.index]
    if len(genes_in_set) < 3: return
    N = len(rnk_series)
    hits_mask = rnk_series.index.isin(genes_in_set)
    rank_metric = rnk_series.values
    Nr = float(np.sum(np.abs(rank_metric[hits_mask])))
    if Nr == 0: return
    running_sum = np.zeros(N)
    miss = (1 - hits_mask).astype(float) / (N - hits_mask.sum())
    hit = hits_mask.astype(float) / Nr * np.abs(rank_metric)
    running_sum = np.cumsum(hit - miss)
    fig, axes = plt.subplots(3, 1, figsize=(7.5, 5.5), sharex=True,
                              gridspec_kw=dict(height_ratios=[3.2, 1.0, 1.2]))
    ax = axes[0]
    ax.plot(running_sum, color="#3B82C4", lw=1.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_ylabel("Running ES")
    ax.set_title(f"{set_name}\nNES={nes:.2g}  nominal p={pval:.3g}  FDR q={fdr:.3g}  size={len(genes_in_set)}")
    ax2 = axes[1]
    hits_idx = np.where(hits_mask)[0]
    ax2.scatter(hits_idx, np.ones_like(hits_idx), s=20, color="#E07A5F", marker="|")
    ax2.set_yticks([]); ax2.set_ylabel("Hits", fontsize=9)
    ax3 = axes[2]
    ax3.bar(range(N), rank_metric, color=["#E07A5F" if h else "#BFBFBF" for h in hits_mask], width=1.0, linewidth=0)
    ax3.set_ylabel("Rank metric", fontsize=9)
    ax3.set_xlabel(f"Ranked genes ({N})")
    fig.suptitle(f"GSEA — {set_name}  (PDO vs PDC)", y=1.0, fontsize=11)
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(GSEA / f"{out_base}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame({"rank": np.arange(N), "gene_symbol": rnk_series.index,
                  "rank_metric": rank_metric, "running_ES": running_sum, "in_set": hits_mask}).to_csv(
        GSEA / f"{out_base}.csv", index=False)

focus_results = []
if not all_gsea.empty:
    for set_name in FOCUS:
        match = all_gsea[all_gsea["Term"] == set_name]
        if match.empty: continue
        row = match.iloc[0]
        gene_set = GENE_SETS[set_name]["genes"]
        nes = float(row["NES"]); pval = float(row["NOM p-val"]); fdr = float(row["FDR q-val"])
        leading = str(row.get("Lead_genes", "")).split(";")
        leading_edge = leading if isinstance(leading, list) else []
        safe = set_name.replace(" ","_").replace("(","").replace(")","").replace("+","_").replace("/","_")
        try: gsea_plot(rnk, gene_set, set_name, f"gsea_{safe}", nes, pval, fdr, leading_edge)
        except Exception as e: print("gsea plot fail:", set_name, e)
        focus_results.append({"pathway": set_name, "NES": nes, "nominal_p": pval,
                              "FDR_q": fdr, "leading_edge_count": len(leading_edge),
                              "leading_edge_genes": ";".join(leading_edge)})

le_df = pd.DataFrame(focus_results)
le_df.to_csv(GSEA / "leading_edge_focused.tsv", sep="\t", index=False)
le_df.to_csv(GSEA / "leading_edge_focused.csv", index=False)

# NES summary
if not le_df.empty:
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    sub = le_df.copy(); sub["dir"] = np.where(sub["NES"]>0, "up_in_PDO", "down_in_PDO")
    colors = ["#E07A5F" if d=="up_in_PDO" else "#3B82C4" for d in sub["dir"]]
    ax.barh(sub["pathway"], sub["NES"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Normalized Enrichment Score (PDO vs PDC)")
    ax.set_title("Focused GSEA — NES (positive = enriched in PDO)")
    for i, r in sub.iterrows():
        ax.text(r["NES"]+(0.05 if r["NES"]>0 else -0.05), i,
                f"p={r['nominal_p']:.2g}, FDR={r['FDR_q']:.2g}",
                va="center", ha="left" if r["NES"]>0 else "right", fontsize=8)
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(GSEA / f"gsea_focused_NES_summary.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    sub.to_csv(GSEA / "gsea_focused_NES_summary.csv", index=False)

# ====== GO / KEGG ORA via Enrichr ======
sig_up = res_df[(res_df["FDR"]<0.05) & (res_df["log2FC"]>1)]["gene_symbol"].dropna().tolist()
sig_down = res_df[(res_df["FDR"]<0.05) & (res_df["log2FC"]<-1)]["gene_symbol"].dropna().tolist()
bg = res_df.dropna(subset=["FDR"])["gene_symbol"].dropna().tolist()

def enrichr_run(gene_list, lib_name, tag):
    if not gene_list: return pd.DataFrame()
    try:
        enr = gp.enrichr(gene_list=gene_list, gene_sets=lib_name, organism="human",
                         outdir=None, no_plot=True)
        r = enr.results
        if r is None or r.empty: return pd.DataFrame()
        r = r.sort_values("Adjusted P-value")
        r.to_csv(GK / f"enrichr_{tag}.tsv", sep="\t", index=False)
        return r
    except Exception as e:
        print(f"enrichr fail {lib_name} {tag}: {e}")
        return pd.DataFrame()

ora_up_bp = enrichr_run(sig_up, "GO_Biological_Process_2023", "up_GO_BP")
ora_up_cc = enrichr_run(sig_up, "GO_Cellular_Component_2023", "up_GO_CC")
ora_up_mf = enrichr_run(sig_up, "GO_Molecular_Function_2023", "up_GO_MF")
ora_up_kegg = enrichr_run(sig_up, "KEGG_2021_Human", "up_KEGG")
ora_down_bp = enrichr_run(sig_down, "GO_Biological_Process_2023", "down_GO_BP")
ora_down_kegg = enrichr_run(sig_down, "KEGG_2021_Human", "down_KEGG")

# Plots: dotplot, lollipop, dumbbell
def plot_ora(df, tag, top=15, plot_type="dot"):
    if df is None or df.empty: return
    d = df.head(top).copy()
    d["-log10_FDR"] = -np.log10(d["Adjusted P-value"].clip(lower=1e-30))
    d["GeneRatio"] = d["Overlap"].apply(lambda x: float(x.split("/")[0])/float(x.split("/")[1]))
    if plot_type == "dot":
        fig, ax = plt.subplots(figsize=(8.5, max(3, 0.42*len(d))))
        sc = ax.scatter(d["GeneRatio"], range(len(d))[::-1], s=d["GeneRatio"]*300+30,
                        c=d["-log10_FDR"], cmap="viridis", edgecolor="black", linewidths=0.4)
        ax.set_yticks(range(len(d))[::-1]); ax.set_yticklabels(d["Term"], fontsize=9)
        ax.set_xlabel("GeneRatio (overlap)"); ax.set_title(f"{tag} — top {len(d)}")
        plt.colorbar(sc, ax=ax, label="-log10(FDR)")
        fig.tight_layout()
        for ext in ("png","pdf","svg"): fig.savefig(GK / f"dotplot_{tag}.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
    elif plot_type == "lollipop":
        fig, ax = plt.subplots(figsize=(8.5, max(3, 0.42*len(d))))
        ax.barh(range(len(d))[::-1], d["GeneRatio"], color="#BFD7EA", edgecolor="grey", height=0.7)
        ax.scatter(d["GeneRatio"], range(len(d))[::-1], s=d["GeneRatio"]*300+30,
                   c=d["-log10_FDR"], cmap="viridis", edgecolor="black", linewidths=0.5, zorder=3)
        ax.set_yticks(range(len(d))[::-1]); ax.set_yticklabels(d["Term"], fontsize=9)
        ax.set_xlabel("GeneRatio"); ax.set_title(f"{tag} — lollipop = -log10(FDR)")
        fig.tight_layout()
        for ext in ("png","pdf","svg"): fig.savefig(GK / f"lollipop_{tag}.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
        d.to_csv(GK / f"lollipop_{tag}.csv", index=False)

for df, tag in [(ora_up_bp, "up_GO_BP"), (ora_up_cc, "up_GO_CC"), (ora_up_mf, "up_GO_MF"),
                (ora_up_kegg, "up_KEGG"), (ora_down_bp, "down_GO_BP"), (ora_down_kegg, "down_KEGG")]:
    plot_ora(df, tag, top=15, plot_type="dot")
    plot_ora(df, tag, top=15, plot_type="lollipop")

# Dumbbell: compare up vs down enrichment for shared terms
def dumbbell_up_vs_down(up_df, down_df, tag, top=10):
    if up_df is None or up_df.empty or down_df is None or down_df.empty: return
    common_terms = set(up_df["Term"]).intersection(down_df["Term"])
    if not common_terms: return
    u = up_df[up_df["Term"].isin(common_terms)].set_index("Term")
    d_ = down_df[down_df["Term"].isin(common_terms)].set_index("Term")
    rows = []
    for t in common_terms:
        u_gr = float(u.loc[t, "Overlap"].split("/")[0])/float(u.loc[t, "Overlap"].split("/")[1])
        d_gr = float(d_.loc[t, "Overlap"].split("/")[0])/float(d_.loc[t, "Overlap"].split("/")[1])
        rows.append({"Term": t, "GeneRatio_up_in_PDO": u_gr,
                     "GeneRatio_down_in_PDO": d_gr,
                     "FDR_up_in_PDO": float(u.loc[t, "Adjusted P-value"]),
                     "FDR_down_in_PDO": float(d_.loc[t, "Adjusted P-value"])})
    df = pd.DataFrame(rows).sort_values("GeneRatio_up_in_PDO", ascending=False).head(top)
    fig, ax = plt.subplots(figsize=(9, max(3, 0.5*len(df))))
    for i, r in df.iterrows():
        y = len(df)-1 - list(df.index).index(i)
        ax.plot([r["GeneRatio_down_in_PDO"], r["GeneRatio_up_in_PDO"]], [y, y], color="grey", lw=2)
        ax.scatter(r["GeneRatio_down_in_PDO"], y, s=80, color="#3B82C4", edgecolor="black", zorder=3, label="down_in_PDO" if i==df.index[0] else "")
        ax.scatter(r["GeneRatio_up_in_PDO"], y, s=80, color="#E07A5F", edgecolor="black", zorder=3, label="up_in_PDO" if i==df.index[0] else "")
    ax.set_yticks(range(len(df))); ax.set_yticklabels(df["Term"], fontsize=9)
    ax.set_xlabel("GeneRatio"); ax.set_title(f"{tag} — dumbbell (up_in_PDO vs down_in_PDO)")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(GK / f"dumbbell_{tag}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    df.to_csv(GK / f"dumbbell_{tag}.csv", index=False)

dumbbell_up_vs_down(ora_up_bp, ora_down_bp, "GO_BP", top=10)
dumbbell_up_vs_down(ora_up_kegg, ora_down_kegg, "KEGG", top=10)

# Combined multi-panel for top GO BP terms across up+down (concatenated)
if not ora_up_bp.empty or not ora_down_bp.empty:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    if not ora_up_bp.empty:
        u = ora_up_bp.head(10).copy()
        u["-log10_FDR"] = -np.log10(u["Adjusted P-value"].clip(lower=1e-30))
        u["GeneRatio"] = u["Overlap"].apply(lambda x: float(x.split("/")[0])/float(x.split("/")[1]))
        axes[0].barh(range(len(u))[::-1], u["GeneRatio"], color="#E07A5F", edgecolor="grey")
        axes[0].set_yticks(range(len(u))[::-1]); axes[0].set_yticklabels(u["Term"], fontsize=9)
        axes[0].set_xlabel("GeneRatio (up_in_PDO)")
        axes[0].set_title("Top GO BP — up_in_PDO")
        for i, r in u.iterrows(): axes[0].text(r["GeneRatio"]+0.005, len(u)-1-i, f"FDR={r['Adjusted P-value']:.2g}", fontsize=7)
    if not ora_down_bp.empty:
        d_ = ora_down_bp.head(10).copy()
        d_["GeneRatio"] = d_["Overlap"].apply(lambda x: float(x.split("/")[0])/float(x.split("/")[1]))
        axes[1].barh(range(len(d_))[::-1], d_["GeneRatio"], color="#3B82C4", edgecolor="grey")
        axes[1].set_yticks(range(len(d_))[::-1]); axes[1].set_yticklabels(d_["Term"], fontsize=9)
        axes[1].set_xlabel("GeneRatio (down_in_PDO / up_in_PDC)")
        axes[1].set_title("Top GO BP — up_in_PDC")
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(GK / f"GO_BP_multi_panel.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)

print("STEP 4 (scores + GSEA + GO/KEGG) complete")