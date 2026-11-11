# Step 3: ZIKV vs Control differential expression of AQP4 / DMD / DAG1 per region x cell type
# Primary test: pseudobulk (sum counts per animal x region x cell type) + edgeR quasi-likelihood,
#   run on the whole transcriptome so FDR is a proper genome-wide FDR.
# Secondary (descriptive only): cell-level Wilcoxon — inflated by pseudoreplication (cells are not independent).
suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
fig_dir <- file.path(proj, "results/figures"); tab_dir <- file.path(proj, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

so <- readRDS(file.path(proj, "data/seurat_working_hvg_annotated.rds"))
genes <- c("AQP4", "DMD", "DAG1")   # DP71 is a DMD isoform — not resolvable from gene-level counts
cond_cols <- c(Control = "#4C78A8", ZIKV = "#E45756")
MIN_CELLS <- 10                      # min cells for a pseudobulk sample to be kept

so <- subset(so, celltype != "Low quality/doublet")
so$celltype <- droplevels(so$celltype); so$celltype_fine <- droplevels(so$celltype_fine)

# ---------------- pseudobulk (full transcriptome, streamed from per-sample checkpoints) ----------------
md <- so@meta.data
md$pb <- paste(md$animal, md$condition, md$region, md$celltype_fine, sep = "|")   # finest level; broad = sum of fine
pb_l <- list()
for (f in list.files(file.path(proj, "data/per_sample_ckpt"), full.names = TRUE)) {
  cnt <- LayerData(readRDS(f)$so, layer = "counts")
  cnt <- cnt[, colnames(cnt) %in% rownames(md), drop = FALSE]
  grp <- factor(md[colnames(cnt), "pb"])
  pb_l[[f]] <- as.matrix(cnt %*% Matrix::sparse.model.matrix(~ 0 + grp))
  colnames(pb_l[[f]]) <- levels(grp)
  rm(cnt); gc(verbose = FALSE)
}
pb_fine <- do.call(cbind, pb_l); rm(pb_l)
info_fine <- md |> count(pb, animal, condition, region, celltype, celltype_fine, name = "n_cells")
info_fine <- info_fine[match(colnames(pb_fine), info_fine$pb), ]
# collapse to broad cell types
info_fine$pb_broad <- paste(info_fine$animal, info_fine$condition, info_fine$region, info_fine$celltype, sep = "|")
agg <- Matrix::sparse.model.matrix(~ 0 + factor(pb_broad), data = info_fine)
pb_broad <- pb_fine %*% as.matrix(agg); colnames(pb_broad) <- sub("^factor\\(pb_broad\\)", "", colnames(agg))
info_broad <- info_fine |> group_by(pb = pb_broad, animal, condition, region, celltype) |>
  summarise(n_cells = sum(n_cells), .groups = "drop")
info_broad <- info_broad[match(colnames(pb_broad), info_broad$pb), ]
saveRDS(list(counts_broad = pb_broad, info_broad = info_broad, counts_fine = pb_fine, info_fine = info_fine),
        file.path(proj, "data/pseudobulk_animal_region_celltype.rds"))
write.csv(info_broad, file.path(tab_dir, "pseudobulk_cell_counts.csv"), row.names = FALSE)

run_pb <- function(pb_mat, pb_info, level, out_prefix) {
  res_all <- list(); cpm_all <- list()
  for (rg in sort(unique(pb_info$region))) for (ct in sort(unique(as.character(pb_info[[level]])))) {
    idx <- which(pb_info$region == rg & pb_info[[level]] == ct & pb_info$n_cells >= MIN_CELLS)
    if (!length(idx)) next
    info <- pb_info[idx, ]
    ok <- sum(info$condition == "Control") >= 2 && sum(info$condition == "ZIKV") >= 2
    y <- DGEList(as.matrix(pb_mat[, idx, drop = FALSE]),
                 group = factor(info$condition, levels = c("Control", "ZIKV")))
    y <- normLibSizes(y)
    cpm_g <- cpm(y, log = FALSE)[genes, , drop = FALSE]
    cpm_all[[paste(rg, ct)]] <- data.frame(info[, c("animal", "condition", "region", "n_cells")], celltype = ct,
                                           t(cpm_g), check.names = FALSE)
    if (!ok) { message("skip ", rg, " / ", ct, " (insufficient replicates)"); next }
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
    design <- model.matrix(~ group, data = y$samples)
    fit <- glmQLFit(estimateDisp(y, design), design, robust = TRUE)
    tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
    tt$gene <- rownames(tt)
    tgt <- data.frame(region = rg, celltype = ct, gene = genes,
                      n_ctrl = sum(info$condition == "Control"), n_zikv = sum(info$condition == "ZIKV"),
                      cells_ctrl = sum(info$n_cells[info$condition == "Control"]),
                      cells_zikv = sum(info$n_cells[info$condition == "ZIKV"]),
                      tested = genes %in% tt$gene)
    tgt <- left_join(tgt, tt[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")], by = "gene")
    res_all[[paste(rg, ct)]] <- tgt
    if (!is.null(out_prefix))
      write.csv(tt[order(tt$PValue), ], file.path(tab_dir, "genome_wide",
                sprintf("%s_%s_%s.csv", out_prefix, rg, gsub("[^A-Za-z0-9]+", "_", ct))), row.names = FALSE)
  }
  res <- bind_rows(res_all)
  # Targeted correction: BH over the 3 genes x all region/celltype tests (hypothesis-driven family)
  res$FDR_targeted <- NA
  res$FDR_targeted[!is.na(res$PValue)] <- p.adjust(res$PValue[!is.na(res$PValue)], "BH")
  list(res = arrange(res, gene, region, celltype), cpm = bind_rows(cpm_all))
}
dir.create(file.path(tab_dir, "genome_wide"), showWarnings = FALSE)
r_broad <- run_pb(pb_broad, info_broad, "celltype", "edgeR_pseudobulk")
r_fine <- run_pb(pb_fine, info_fine, "celltype_fine", NULL)
res <- r_broad$res; cpm_df <- r_broad$cpm
write.csv(res, file.path(tab_dir, "target_genes_pseudobulk_edgeR_ZIKV_vs_Control.csv"), row.names = FALSE)
write.csv(cpm_df, file.path(tab_dir, "target_genes_pseudobulk_CPM.csv"), row.names = FALSE)
write.csv(r_fine$res, file.path(tab_dir, "target_genes_pseudobulk_edgeR_ZIKV_vs_Control_FINE_subtypes.csv"), row.names = FALSE)

# ---------------- secondary: cell-level Wilcoxon ----------------
so$grp <- paste(so$region, so$celltype, sep = "|")
wil <- list()
for (g in unique(so$grp)) {
  sub <- subset(so, grp == g)
  tb <- table(sub$condition)
  if (length(tb) < 2 || any(tb < 20)) next
  fm <- FindMarkers(sub, ident.1 = "ZIKV", ident.2 = "Control", group.by = "condition",
                    features = genes, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
  fm$gene <- rownames(fm); fm$region <- sub$region[1]; fm$celltype <- sub$celltype[1]
  wil[[g]] <- fm
}
wil <- bind_rows(wil) |> select(region, celltype, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) |>
  rename(pct_ZIKV = pct.1, pct_Control = pct.2) |> arrange(gene, region, celltype)
write.csv(wil, file.path(tab_dir, "target_genes_cell_level_wilcoxon_DESCRIPTIVE_ONLY.csv"), row.names = FALSE)

# store results in the object
so@misc$target_gene_DE <- list(pseudobulk_edgeR = res, pseudobulk_edgeR_fine = r_fine$res,
                                cell_level_wilcoxon = wil, pseudobulk_cpm = cpm_df)
so$grp <- NULL
saveRDS(so, file.path(proj, "data/seurat_working_hvg_annotated.rds"))
print(res |> filter(!is.na(PValue)) |> arrange(PValue) |> head(20))
message("done — figures: scripts/03b_figures.R")
