"""Step 5f — GSEA + GO/KEGG for Ulixertinib vs Ctrl.

Reuses the DE result from 05e (vendor counts). Runs prerank GSEA on
Hallmark + curated sets, then ORA via Enrichr on up/down DE gene lists.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp

TURN = Path("/home/susu/Rprojects/PCP202607/outputs/turn_20260802065502_636832")
A2 = TURN / "03_Ulixertinib_vs_Ctrl"
GE = A2 / "GSEA"; GE.mkdir(exist_ok=True)
GK = A2 / "GO_KEGG"; GK.mkdir(exist_ok=True)

CURATED = {
    "Stemness (Wong 2008 + Bhatt 2012)": [
        "POU5F1","SOX2","NANOG","MYC","KLF4","BMI1","CD44","ALDH1A1","PROM1","NES","ITGA6","NOTCH1",
        "WNT5A","FZD2","FGF2","FGF5","GATA6","LGR5","LGR4","ABCG2","CDH1","VIM","SNAI1","SNAI2",
        "ZEB1","ZEB2","TWIST1","CDH2","CHI3L1","FN1","LAMA5","ITGB1","CD9","CD63","CD81","TSG101",
    ],
    "Proliferation (Hallmark E2F + G2M)": [
        "MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCNA2","CCNE1","CCNE2","CCND1","CCND2","CCND3",
        "CDK1","CDK2","CDK4","CDK6","E2F1","E2F2","E2F3","E2F4","E2F5","E2F6","E2F7","E2F8",
        "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM10","CDC6","ORC1","ORC2","ORC3","ORC4","ORC5","ORC6",
        "CDC20","BUB1","BUB1B","MAD2L1","BIRC5","AURKA","AURKB","PLK1","TTK","KIF2C","KIF11","KIF14","KIF20A",
        "CENPA","CENPB","CENPE","CENPF","NDC80","NUF2","SPC24","SPC25","RAD21","SMC1A","SMC3","STAG1","STAG2",
        "RRM1","RRM2","TYMS","DHFR","POLA1","POLD1","POLD2","POLE","POLE2","FEN1","LIG1","UBE2C","UBE2S",
        "TUBB","TUBA1B","TUBB4B","H2AX","H2AZ1","FOXM1","MYBL2","BRCA1","BRCA2","RAD51","CHEK1","CHEK2",
        "TP53","MDM2","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","RB1","RBL1","RBL2",
    ],
    "Apoptosis (Hallmark Apoptosis core)": [
        "CASP3","CASP6","CASP7","CASP8","CASP9","CASP10","CASP4","CASP1","CASP2","CASP5","CASP14",
        "FAS","FASLG","FADD","TRADD","TNFRSF10A","TNFRSF10B","TNFRSF1A","TNFRSF1B","TNFA","TNF",
        "BID","BAX","BAK1","BCL2","BCL2L1","BCL2L11","BCL2L2","MCL1","BIRC2","BIRC3","BIRC5","XIAP",
        "CYCS","APAF1","AIFM1","ENDOG","PARP1","PARP2","DFFA","DFFB","TP53","TP63","TP73",
        "TRAF2","TRAF3","TRAF6","RIPK1","RIPK3","MLKL","JAK1","JAK2","STAT1","STAT3","NFKB1","NFKBIA",
    ],
    "MAPK/ERK cascade (KEGG hsa04010 + Reactome RAF/MAPK)": [
        "MAPK1","MAPK3","MAP2K1","MAP2K2","MAP3K1","MAP3K2","MAP3K3","MAP3K4","MAP3K5","MAP3K11","MAP3K12",
        "RAF1","ARAF","BRAF","KRAS","NRAS","HRAS","SOS1","SOS2","GRB2","PTPN11","RASGRP1",
        "ELK1","ETS1","ETS2","JUN","FOS","ATF2","MEF2C","MSK1","MSK2","MNK1","MNK2","RSK1","RSK2",
        "PLA2G4A","PLCB1","PLCB2","PLCD1","PLCG1","PLCG2","PRKCA","PRKCB","PRKCG","PRKD1","PRKD2","PRKD3",
        "PDGFRA","PDGFRB","PDGFA","PDGFB","VEGFA","VEGFB","VEGFC","FLT1","KDR","FLT4","EGFR","ERBB2",
        "FGFR1","FGFR2","FGFR3","NGFR","TRKA","NTRK1","NTRK2","IGF1R","INSR","MET","KIT",
    ],
    "Cell cycle (GO 0007049 + Reactome Cell Cycle)": [
        "CDK1","CDK2","CDK4","CDK6","CDK7","CDK9","CCNB1","CCNB2","CCNA2","CCNE1","CCNE2","CCND1","CCND2","CCND3",
        "E2F1","E2F2","E2F3","E2F4","E2F5","E2F6","E2F7","E2F8","RB1","RBL1","RBL2","TP53","MDM2","MDM4",
        "CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDC25A","CDC25B","CDC25C","CHEK1","CHEK2","ATR","ATM",
        "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM10","PCNA","ORC1","ORC2","ORC3","ORC4","ORC5","ORC6","CDC6",
        "FOXM1","MYBL2","BIRC5","AURKA","AURKB","PLK1","TTK","CENPA","CENPB","CENPE","NDC80","SPC24","MAD2L1","BUB1","BUB1B",
        "BRCA1","BRCA2","RAD51","RFC1","RFC2","RFC3","RFC4","RFC5","LIG1","FEN1","TYMS","RRM1","RRM2",
    ],
}

# Load DE results (from 05e vendor pipeline)
res = pd.read_csv(A2 / "differential_expression" / "all_genes_Ulixertinib_vs_Ctrl.tsv", sep="\t")

rank_df = res.dropna(subset=["stat"]).copy()
rank_df = rank_df.sort_values("stat", ascending=False)
ranked_out = rank_df[["gene_symbol","stat"]].rename(columns={"gene_symbol":"#Name","stat":"score"})
ranked_out.to_csv(GE / "prerank_input.tsv", sep="\t", index=False)
rnk = ranked_out.set_index("#Name")["score"]

# Curated
res_curated = gp.prerank(rnk=rnk, gene_sets=CURATED,
                         permutation_num=1000, min_size=5, max_size=2000,
                         seed=42, no_plot=True, verbose=False)
curated_df = res_curated.res2d.copy()
curated_df["library"] = "Curated (5 focused)"
curated_df.to_csv(GE / "prerank_curated.tsv", sep="\t", index=False)
print(f"Curated prerank: {len(curated_df)} pathways")

# Hallmark
hallmark = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
res_hall = gp.prerank(rnk=rnk, gene_sets=hallmark, permutation_num=1000,
                      min_size=10, max_size=2000, seed=42, no_plot=True, verbose=False)
hall_df = res_hall.res2d.copy()
hall_df["library"] = "MSigDB_Hallmark_2020"
hall_df.to_csv(GE / "prerank_hallmark.tsv", sep="\t", index=False)
print(f"Hallmark prerank: {len(hall_df)} pathways")

# Combined
all_gsea = pd.concat([curated_df, hall_df], ignore_index=True)
all_gsea.to_csv(GE / "gsea_all_pathways.tsv", sep="\t", index=False)

# NES summary for focused
focused = curated_df.copy()
focused_summary = pd.DataFrame({
    "pathway": focused["Term"],
    "NES": focused["NES"],
    "nominal_p": focused["NOM p-val"],
    "FDR_q": focused["FDR q-val"],
    "leading_edge_count": focused["Lead_genes"].fillna("").str.split(";").apply(lambda x: len([g for g in x if g])),
    "leading_edge_genes": focused["Lead_genes"],
    "dir": focused["NES"].apply(lambda v: "up_in_Ulixertinib" if v>0 else "down_in_Ulixertinib"),
})
focused_summary.to_csv(GE / "gsea_focused_NES_summary.csv", index=False)

le_rows = []
for _, r in focused.iterrows():
    le = [g for g in str(r["Lead_genes"]).split(";") if g]
    le_rows.append({"pathway": r["Term"], "NES": r["NES"], "leading_edge_gene": ";".join(le), "n_leading_edge": len(le)})
le_df = pd.DataFrame(le_rows)
le_df.to_csv(GE / "leading_edge_focused.tsv", sep="\t", index=False)

# NES bar
fig, ax = plt.subplots(figsize=(7, 4))
sorted_f = focused.sort_values("NES")
colors = ["#4DBBD5" if x>0 else "#E64B35" for x in sorted_f["NES"]]
ax.barh(sorted_f["Term"], sorted_f["NES"], color=colors, edgecolor="k")
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("NES (positive = enriched in Ulixertinib)")
ax.set_title("GSEA NES summary — Ulixertinib vs Ctrl (5 focused pathways)")
fig.tight_layout()
for ext in ("png","pdf","svg"): fig.savefig(GE / f"gsea_focused_NES_summary.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

# Single-pathway GSEA plots for the 5 focused
for name in CURATED.keys():
    safe = "".join(c if c.isalnum() or c=="_" else "_" for c in name)
    sets = {name: CURATED[name]}
    try:
        gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=5, max_size=2000,
                   seed=42, outdir=str(GE / f"_prerank_{safe}"), no_plot=False, verbose=False)
        # gseapy writes to outdir/<safename>.png and <safename>.enplot.png; rename to canonical
        outd = GE / f"_prerank_{safe}"
        if (outd / f"{safe}.png").exists():
            import shutil
            shutil.copy(outd / f"{safe}.png", GE / f"gsea_{safe}.png")
        if (outd / f"{safe}.enplot.png").exists():
            import shutil
            shutil.copy(outd / f"{safe}.enplot.png", GE / f"gsea_{safe}_running.png")
    except Exception as e:
        print(f"GSEA plot for {name} failed: {e}")
print("GSEA done")

# === GO/KEGG ORA ===
sig_up = res[res["direction"]=="up_in_Ulixertinib"]["gene_symbol"].dropna().tolist()
sig_dn = res[res["direction"]=="up_in_Ctrl"]["gene_symbol"].dropna().tolist()
print(f"up_in_Ulixertinib genes: {len(sig_up)}; up_in_Ctrl genes: {len(sig_dn)}")

def run_enrichr(gene_list, tag):
    if len(gene_list) < 5:
        return None
    libs = {"GO_BP":"GO_Biological_Process_2023",
            "GO_CC":"GO_Cellular_Component_2023",
            "GO_MF":"GO_Molecular_Function_2023",
            "KEGG":"KEGG_2021_Human"}
    frames = []
    for short, full in libs.items():
        try:
            enr = gp.enrichr(gene_list=gene_list, gene_sets=full, organism="human", outdir=None, no_plot=True)
            edf = enr.results.copy()
            edf = edf.sort_values("Adjusted P-value")
            edf["library_short"] = short
            edf.to_csv(GK / f"enrichr_{tag}_{short}.tsv", sep="\t", index=False)
            frames.append(edf)
        except Exception as e:
            print(f"Enrichr {short} failed: {e}")
            frames.append(pd.DataFrame())
    return pd.concat(frames, ignore_index=True) if frames else None

up_enr = run_enrichr(sig_up, "up")
dn_enr = run_enrichr(sig_dn, "down")

def render_go_kegg(enr_df, contrast, palette):
    if enr_df is None or enr_df.empty: return
    for short in enr_df["library_short"].unique():
        sub = enr_df[enr_df["library_short"] == short].head(10).copy()
        sub["-log10_FDR"] = -np.log10(sub["Adjusted P-value"].clip(lower=1e-30))
        # Dot plot
        fig, ax = plt.subplots(figsize=(7, 4.5))
        sizes = (sub["Overlap"].str.split("/").str[0].astype(int) /
                 sub["Overlap"].str.split("/").str[1].astype(int) * 200 + 20)
        ax.scatter(sub["-log10_FDR"], range(len(sub))[::-1], s=sizes, color=palette, edgecolor="k")
        ax.set_yticks(range(len(sub))[::-1])
        ax.set_yticklabels(sub["Term"], fontsize=9)
        ax.set_xlabel("-log10(FDR)")
        ax.set_title(f"{contrast} {short} top 10 — Ulixertinib vs Ctrl")
        fig.tight_layout()
        for ext in ("png","pdf","svg"): fig.savefig(GK / f"dotplot_{contrast}_{short}.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
        sub.to_csv(GK / f"dotplot_{contrast}_{short}.csv", index=False)
        # Lollipop
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.hlines(y=range(len(sub))[::-1], xmin=0, xmax=sub["-log10_FDR"], color=palette, lw=2)
        ax.scatter(sub["-log10_FDR"], range(len(sub))[::-1], s=60, color=palette, edgecolor="k")
        ax.set_yticks(range(len(sub))[::-1])
        ax.set_yticklabels(sub["Term"], fontsize=9)
        ax.set_xlabel("-log10(FDR)")
        ax.set_title(f"{short} {contrast} — Ulixertinib vs Ctrl")
        fig.tight_layout()
        for ext in ("png","pdf","svg"): fig.savefig(GK / f"lollipop_{contrast}_{short}.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
        sub.to_csv(GK / f"lollipop_{contrast}_{short}.csv", index=False)

render_go_kegg(up_enr, "up", "#D6604D")
render_go_kegg(dn_enr, "down", "#4393C3")

# Dumbbell
if up_enr is not None and dn_enr is not None and not up_enr.empty and not dn_enr.empty:
    for short in ["GO_BP","KEGG"]:
        u = up_enr[up_enr["library_short"]==short][["Term","Adjusted P-value"]].rename(columns={"Adjusted P-value":"up_FDR"})
        d = dn_enr[dn_enr["library_short"]==short][["Term","Adjusted P-value"]].rename(columns={"Adjusted P-value":"down_FDR"})
        merged = u.merge(d, on="Term", how="inner").head(15)
        if merged.empty: continue
        merged["up_-log10"] = -np.log10(merged["up_FDR"].clip(lower=1e-30))
        merged["down_-log10"] = -np.log10(merged["down_FDR"].clip(lower=1e-30))
        fig, ax = plt.subplots(figsize=(8, 5))
        y = range(len(merged))[::-1]
        ax.hlines(y, merged["down_-log10"], merged["up_-log10"], color="grey", lw=1)
        ax.scatter(merged["up_-log10"], y, color="#D6604D", s=60, edgecolor="k", label="Up in Ulixertinib")
        ax.scatter(merged["down_-log10"], y, color="#4393C3", s=60, edgecolor="k", label="Up in Ctrl (down in Ulixertinib)")
        ax.set_yticks(y); ax.set_yticklabels(merged["Term"], fontsize=8)
        ax.set_xlabel("-log10(FDR)")
        ax.set_title(f"{short} up-vs-down — Ulixertinib vs Ctrl")
        ax.legend()
        fig.tight_layout()
        for ext in ("png","pdf","svg"): fig.savefig(GK / f"dumbbell_{short}.{ext}", dpi=600, bbox_inches="tight")
        plt.close(fig)
        merged.to_csv(GK / f"dumbbell_{short}.csv", index=False)

print("GO/KEGG done")
print("DONE step 5f")