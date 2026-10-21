args <- commandArgs(trailingOnly = TRUE)
contrast <- args[1]  # "PDO_vs_PDC" or "Ulixertinib_vs_Ctrl"
work <- args[2]
out_dir <- file.path(work, "09_GSEA_visualization")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

cat("=== GSEA contrast:", contrast, "===\n")

# ----------------------------------------------------------------
# 1. 构造 symbol -> score (rank) 向量
# ----------------------------------------------------------------
# PDC/PDO: 必须从 DE 表里取 (log2FC)
# Ulixertinib: 已有 prerank_input.tsv (symbol, score)
if (contrast == "PDO_vs_PDC") {
  de_path <- file.path(work, "02_PDC_vs_PDO", "differential_expression",
                       "all_genes_PDO_vs_PDC.tsv")
  de <- read.delim(de_path, stringsAsFactors = FALSE)
  stopifnot("gene_symbol" %in% colnames(de))
  de <- de[!is.na(de$gene_symbol) & de$gene_symbol != "" & !is.na(de$log2FC), ]
  de <- de[!duplicated(de$gene_symbol), ]
  genelist <- de$log2FC
  names(genelist) <- de$gene_symbol
  genelist <- sort(genelist, decreasing = TRUE)
} else {
  pr_path <- file.path(work, "03_Ulixertinib_vs_Ctrl", "GSEA", "prerank_input.tsv")
  pr <- read.delim(pr_path, stringsAsFactors = FALSE, comment.char = "#")
  if (!"Name" %in% colnames(pr)) colnames(pr)[1] <- "Name"
  if (!"score" %in% colnames(pr)) colnames(pr)[2] <- "score"
  pr <- pr[, c("Name","score")]
  pr <- pr[!is.na(pr$score) & !is.na(pr$Name) & pr$Name != "", ]
  pr <- pr[!duplicated(pr$Name), ]
  genelist <- pr$score
  names(genelist) <- pr$Name
  genelist <- sort(genelist, decreasing = TRUE)
}
cat("ranked gene list length:", length(genelist), "\n")
cat("first 5:", paste(head(names(genelist),5), head(genelist,5), sep="=", collapse="; "), "\n")

# ----------------------------------------------------------------
# 1b. symbol -> entrez (gseGO 要求 entrez)
# ----------------------------------------------------------------
map_df <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = unique(names(genelist)),
                                keytype = "SYMBOL",
                                columns = c("ENTREZID"))
map_df <- map_df[!is.na(map_df$ENTREZID), ]
map_df <- map_df[!duplicated(map_df$SYMBOL) & !duplicated(map_df$ENTREZID), ]
cat("symbol->entrez mapped:", nrow(map_df), "/", length(genelist), "\n")
genelist_entrez <- genelist[map_df$SYMBOL]
names(genelist_entrez) <- as.character(map_df$ENTREZID)
genelist_entrez <- genelist_entrez[!is.na(names(genelist_entrez))]
genelist_entrez <- sort(genelist_entrez, decreasing = TRUE)
cat("entrez ranked length:", length(genelist_entrez), "\n")

# ----------------------------------------------------------------
# 2. 库: GO:BP + KEGG + Hallmark
# ----------------------------------------------------------------
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- msig_h %>% dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame() %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)
msig_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_term_sets <- msig_kegg %>% dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame() %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

# ----------------------------------------------------------------
# 3. 运行 GSEA
# ----------------------------------------------------------------
cat("\n--- gseGO BP ---\n")
suppressWarnings(suppressMessages({
  ego <- gseGO(geneList = genelist_entrez,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               minGSSize = 15,
               maxGSSize = 500,
               pvalueCutoff = 1,            # 不过滤,画棒棒糖时统一用 FDR
               verbose = FALSE,
               pAdjustMethod = "BH")
}))
cat("  ego @result rows:", nrow(as.data.frame(ego)), "\n")

cat("\n--- gseKEGG ---\n")
suppressWarnings(suppressMessages({
  ekegg <- gseKEGG(geneList = genelist_entrez,
                   organism = "hsa",
                   minGSSize = 15,
                   maxGSSize = 500,
                   pvalueCutoff = 1,
                   verbose = FALSE,
                   pAdjustMethod = "BH",
                   use_internal_data = FALSE)
}))
cat("  ekegg @result rows:", nrow(as.data.frame(ekegg)), "\n")

cat("\n--- fgsea Hallmark ---\n")
# Hallmark 走 fgsea (用 msigdbr Hallmark, 用 symbol genelist)
hallmark_list <- hallmark %>% dplyr::group_by(term) %>%
  dplyr::summarise(genes = list(gene), .groups = "drop") %>%
  dplyr::mutate(genes = lapply(genes, unique)) %>%
  tibble::deframe()
library(fgsea)
efg <- fgsea(pathways = hallmark_list,
             stats = genelist,
             eps = 0,
             minSize = 15,
             maxSize = 500)
efg <- efg %>% as.data.frame()
cat("  fgsea hallmark rows:", nrow(efg), "\n")

# ----------------------------------------------------------------
# 4. 汇总为 1 张大表 (GO + KEGG + Hallmark 全部)
# ----------------------------------------------------------------
go_df <- as.data.frame(ego) %>%
  dplyr::mutate(library = "GO:BP",
                term_id = ID,
                term_name = Description,
                leading_edge = core_enrichment,
                gene_ratio = paste0(setSize, " genes"),
                padj = p.adjust,
                ES = enrichmentScore) %>%
  as.data.frame()
kegg_df <- as.data.frame(ekegg) %>%
  dplyr::mutate(library = "KEGG",
                term_id = ID,
                term_name = Description,
                leading_edge = core_enrichment,
                gene_ratio = paste0(setSize, " genes"),
                padj = p.adjust,
                ES = enrichmentScore) %>%
  as.data.frame()
hallmark_df <- efg %>%
  dplyr::mutate(library = "Hallmark",
                term_id = pathway,
                term_name = pathway,
                leading_edge = leadingEdge,
                gene_ratio = paste0(size, " genes"),
                NES = NES,
                pvalue = pval,
                padj = padj,
                ES = NA) %>%
  dplyr::select(library, term_id, term_name, ES, NES, pvalue, padj,
                gene_ratio, leading_edge)

keep_cols <- c("library","term_id","term_name","ES","NES","pvalue","padj",
               "gene_ratio","leading_edge")
cat("go_df cols:", paste(colnames(go_df), collapse=","), "\n")
cat("kegg_df cols:", paste(colnames(kegg_df), collapse=","), "\n")
cat("hallmark_df cols:", paste(colnames(hallmark_df), collapse=","), "\n")
# 强制所有 df 共有的列 & ES 转 numeric
go_df[,        keep_cols]    <- lapply(go_df[, keep_cols],       as.character)
kegg_df[,      keep_cols]    <- lapply(kegg_df[, keep_cols],     as.character)
hallmark_df[,  keep_cols]    <- lapply(hallmark_df[, keep_cols], as.character)
go_df$ES        <- as.numeric(go_df$ES)
kegg_df$ES      <- as.numeric(kegg_df$ES)
hallmark_df$ES  <- as.numeric(hallmark_df$ES)
go_df$NES       <- as.numeric(go_df$NES)
kegg_df$NES     <- as.numeric(kegg_df$NES)
hallmark_df$NES <- as.numeric(hallmark_df$NES)
go_df$pvalue    <- as.numeric(go_df$pvalue)
kegg_df$pvalue  <- as.numeric(kegg_df$pvalue)
hallmark_df$pvalue <- as.numeric(hallmark_df$pvalue)
go_df$padj      <- as.numeric(go_df$padj)
kegg_df$padj    <- as.numeric(kegg_df$padj)
hallmark_df$padj <- as.numeric(hallmark_df$padj)
all_df <- rbind(go_df[, keep_cols], kegg_df[, keep_cols], hallmark_df[, keep_cols])
all_df <- all_df %>% dplyr::arrange(padj, dplyr::desc(abs(NES)))
all_df$contrast <- contrast
write.table(all_df, file.path(out_dir, paste0("gsea_ALL_", contrast, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("ALL summary rows:", nrow(all_df), " saved\n")

# 各子库单独
write.table(go_df[, keep_cols],    file.path(out_dir, paste0("gsea_GO_", contrast, ".tsv")),    sep="\t", quote=FALSE, row.names=FALSE)
write.table(kegg_df[, keep_cols],  file.path(out_dir, paste0("gsea_KEGG_", contrast, ".tsv")),  sep="\t", quote=FALSE, row.names=FALSE)
write.table(hallmark_df[, keep_cols], file.path(out_dir, paste0("gsea_Hallmark_", contrast, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)

# 存 .RData 供再画
saveRDS(list(ego = ego, ekegg = ekegg, efg = efg, all = all_df),
        file.path(out_dir, paste0("gsea_objects_", contrast, ".RDS")))

cat("\nDONE contrast:", contrast, "\n")
