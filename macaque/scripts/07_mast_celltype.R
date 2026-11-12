# Step 7: cell-level differential expression with MAST, ZIKV vs Control, within each cell type.
#
# This complements the pseudobulk edgeR of 03/05; it does not replace it. MAST models what pseudobulk
# throws away -- the zero-inflated per-cell distribution and the cellular detection rate -- but cells
# from the same animal are not independent replicates, and there are only 3 vs 2 animals here. Two models
# are therefore fitted:
#   (a) hurdle, fixed effects only:  ~ condition + cngeneson + region
#       The conventional MAST call. Treats 176,570 cells as replicates, so its P-values are
#       anticonservative by orders of magnitude. Reported for effect size and detection rate, NOT evidence.
#   (b) hurdle with the animal as a random effect:  ~ condition + cngeneson + region + (1|animal)
#       The defensible test: the animal is the unit of replication. Fitted for the endfoot-anchor genes.
# Both are run on panel genes only, read straight from the per-sample checkpoints, so the 2.5 GB object
# is never loaded. Region is a fixed covariate here (the three regions are pooled), which is why effect
# sizes are not expected to match the per-region pseudobulk exactly.
suppressPackageStartupMessages({
  library(Matrix); library(MAST); library(data.table); library(dplyr); library(SummarizedExperiment)
})
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
tab_dir <- file.path(proj, "results/tables")
set.seed(1)

ANCHOR <- c("DAG1", "DMD", "UTRN", "LARGE1", "SNTA1", "DTNA")          # the endfoot anchor asked for
panel <- list(
  "Endfoot / DGC (astrocyte side)" = c("AQP4", "KCNJ10", "DMD", "DAG1", "SNTA1", "SNTB1", "DTNA", "DTNB",
                                       "UTRN", "MLC1", "HEPACAM", "LARGE1", "POMT1"),
  "Basement-membrane ligands" = c("AGRN", "LAMA2", "LAMA4", "LAMB1", "LAMB2", "LAMC1", "HSPG2", "NID1",
                                  "COL4A1", "COL4A2"),
  "Vascular tone and BBB" = c("AVPR1A", "ACTA2", "MYH11", "PDGFRB", "CLDN5", "OCLN", "TJP1"),
  "Astrocyte reactivity" = c("GFAP", "VIM", "CD44"),
  "Interferon response" = c("ISG15", "IFI6", "IFI27", "MX1", "IFIT1", "IFIT3", "STAT1", "OAS1"))
panel_df <- data.frame(group = rep(names(panel), lengths(panel)), gene = unlist(panel), row.names = NULL)
genes <- unique(panel_df$gene)

MIN_CELLS <- 50      # per condition, within a cell type
MIN_ANIMAL_CELLS <- 10

# ---- panel-only cell-level counts, assembled from the per-sample checkpoints ----
cache <- file.path(proj, "data/panel_cell_expr.rds")
if (!file.exists(cache)) {
  md <- readRDS(file.path(proj, "data/cell_metadata.rds"))
  files <- sort(list.files(file.path(proj, "data/per_sample_ckpt"), full.names = TRUE))
  mats <- lapply(files, function(f) {
    m <- SeuratObject::LayerData(readRDS(f)$so, layer = "counts")
    m <- m[genes, colnames(m) %in% md$cell, drop = FALSE]
    gc(verbose = FALSE); m
  })
  counts <- do.call(cbind, mats); rm(mats); gc()
  md <- md[match(colnames(counts), md$cell), ]
  stopifnot(identical(md$cell, colnames(counts)))
  saveRDS(list(counts = counts, md = md), cache)
  rm(counts, md); gc()
}
pe <- readRDS(cache); counts <- pe$counts; md <- pe$md; rm(pe); gc()
message("panel matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

# MAST expects log2 expression; counts per 10,000 keeps the usual 10x scaling (Seurat's scale.factor)
lognorm <- log2(t(t(counts) / md$nCount_RNA) * 1e4 + 1)
md$cngeneson <- as.numeric(scale(md$nFeature_RNA))                      # cellular detection rate
md$condition <- factor(md$condition, levels = c("Control", "ZIKV"))
md$animal <- factor(md$animal); md$region <- factor(md$region)

# ---- descriptive: where the anchor genes are detected at all, by cell type ----
det <- bind_rows(lapply(split(seq_len(ncol(counts)), md$celltype), function(j) {
  data.frame(celltype = as.character(md$celltype[j[1]]), n_cells = length(j), gene = ANCHOR,
             pct_detected = 100 * rowMeans(counts[ANCHOR, j, drop = FALSE] > 0),
             mean_log2 = rowMeans(lognorm[ANCHOR, j, drop = FALSE]), row.names = NULL)
})) |> arrange(gene, desc(pct_detected))
write.csv(det, file.path(tab_dir, "mast_endfoot_detection_by_celltype.csv"), row.names = FALSE)

# ---- helper: the cell types testable for a given grouping column ----
testable <- function(col) {
  tab <- md |> count(.data[[col]], condition, animal) |> filter(n >= MIN_ANIMAL_CELLS)
  keep <- tab |> count(.data[[col]], condition) |> tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) |>
    filter(Control >= 2, ZIKV >= 2)
  tot <- md |> count(.data[[col]], condition) |> tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) |>
    filter(Control >= MIN_CELLS, ZIKV >= MIN_CELLS)
  intersect(as.character(keep[[col]]), as.character(tot[[col]]))
}

# ---- fit one cell type ----
fit_one <- function(idx, use_genes, random_effect) {
  sub <- md[idx, ]
  g <- use_genes[rowSums(counts[use_genes, idx, drop = FALSE] > 0) >= 10]   # needs some signal to fit
  if (!length(g)) return(NULL)
  sca <- FromMatrix(as.matrix(lognorm[g, idx, drop = FALSE]),
                    cData = data.frame(wellKey = sub$cell, condition = sub$condition, cngeneson = sub$cngeneson,
                                       region = sub$region, animal = sub$animal),
                    fData = data.frame(primerid = g), check_sanity = FALSE)
  form <- if (random_effect) ~ condition + cngeneson + region + (1 | animal) else ~ condition + cngeneson + region
  zz <- try(suppressMessages(suppressWarnings(
    # nAGQ applies to the binomial (discrete) half only; lmer, which fits the continuous half, rejects it
    if (random_effect) zlm(form, sca, method = "glmer", ebayes = FALSE, strictConvergence = FALSE,
                           fitArgsD = list(nAGQ = 0))
    else zlm(form, sca))), silent = TRUE)
  if (inherits(zz, "try-error")) return(NULL)
  lr <- try(suppressMessages(suppressWarnings(lrTest(zz, CoefficientHypothesis("conditionZIKV")))), silent = TRUE)
  if (inherits(lr, "try-error")) return(NULL)
  cf <- try(suppressWarnings(summary(zz, logFC = TRUE)$datatable), silent = TRUE)
  lfc <- rep(NA_real_, length(g)); names(lfc) <- g
  if (!inherits(cf, "try-error")) {
    d <- as.data.table(cf)[contrast == "conditionZIKV" & component == "logFC"]
    if (nrow(d)) lfc[as.character(d$primerid)] <- d$coef
  }
  data.frame(gene = g, logFC_MAST = unname(lfc[g]),
             PValue = lr[g, "hurdle", "Pr(>Chisq)"],
             pct_zikv = 100 * rowMeans(counts[g, idx[sub$condition == "ZIKV"], drop = FALSE] > 0),
             pct_ctrl = 100 * rowMeans(counts[g, idx[sub$condition == "Control"], drop = FALSE] > 0),
             n_cells = length(idx), row.names = NULL)
}

run_model <- function(col, use_genes, random_effect, label) {
  cts <- testable(col)
  message("[", label, "] ", length(cts), " cell types x ", length(use_genes), " genes")
  out <- bind_rows(lapply(cts, function(ct) {
    idx <- which(as.character(md[[col]]) == ct)
    t0 <- Sys.time(); r <- fit_one(idx, use_genes, random_effect)
    message("  ", ct, " (", length(idx), " cells) ",
            if (is.null(r)) "FAILED" else sprintf("%.0fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(r)) NULL else cbind(celltype = ct, r)
  }))
  out$FDR <- p.adjust(out$PValue, "BH")
  out
}

# ---- (a) conventional MAST, full panel, every cell type: effect sizes, not evidence ----
fixed <- run_model("celltype", genes, FALSE, "fixed effects")
fixed <- left_join(fixed, panel_df, by = "gene") |>
  select(celltype, group, gene, logFC_MAST, PValue, FDR, pct_ctrl, pct_zikv, n_cells) |>
  arrange(PValue)
write.csv(fixed, file.path(tab_dir, "mast_panel_ZIKV_vs_Control_by_celltype_FIXED_pseudoreplicated.csv"), row.names = FALSE)

# ---- (b) animal as random effect, endfoot anchor genes: the defensible test ----
mixed <- run_model("celltype", ANCHOR, TRUE, "animal random effect") |>
  select(celltype, gene, logFC_MAST, PValue, FDR, pct_ctrl, pct_zikv, n_cells) |> arrange(PValue)
write.csv(mixed, file.path(tab_dir, "mast_endfoot_anchor_animalRE_by_celltype.csv"), row.names = FALSE)

# astrocyte subtypes too, same model
mixed_fine <- run_model("celltype_fine", ANCHOR, TRUE, "animal RE, fine") |>
  filter(grepl("^Astrocytes", celltype)) |>
  select(celltype_fine = celltype, gene, logFC_MAST, PValue, FDR, pct_ctrl, pct_zikv, n_cells) |> arrange(PValue)
write.csv(mixed_fine, file.path(tab_dir, "mast_endfoot_anchor_animalRE_astrocyte_subtypes.csv"), row.names = FALSE)

# ---- how much does ignoring the animal inflate significance? ----
pb <- read.csv(file.path(tab_dir, "glymphatic_panel_pseudobulk_edgeR_broad.csv")) |>
  filter(gene %in% ANCHOR, tested) |> group_by(celltype, gene) |>
  summarise(edgeR_logFC_median = median(logFC), edgeR_P_min = min(PValue), .groups = "drop")
cmp <- fixed |> filter(gene %in% ANCHOR) |> select(celltype, gene, logFC_MAST_fixed = logFC_MAST, P_fixed = PValue) |>
  full_join(mixed |> select(celltype, gene, logFC_MAST_RE = logFC_MAST, P_animalRE = PValue), by = c("celltype", "gene")) |>
  full_join(pb, by = c("celltype", "gene")) |>
  mutate(P_inflation_log10 = log10(P_animalRE) - log10(P_fixed)) |> arrange(P_animalRE)
write.csv(cmp, file.path(tab_dir, "mast_vs_pseudobulk_endfoot_comparison.csv"), row.names = FALSE)

show <- function(x, n = 40) base::print(tibble::as_tibble(x), n = n, width = 200)
message("\n--- (b) animal random effect, P < 0.05 ---"); show(mixed |> filter(PValue < 0.05))
message("\n--- astrocyte subtypes ---"); show(mixed_fine |> filter(PValue < 0.05))
message("\n--- fixed vs random effect, astrocytes ---"); show(cmp |> filter(celltype == "Astrocytes"))
message("\n--- median P inflation from ignoring the animal (log10 units) ---")
message(round(median(cmp$P_inflation_log10, na.rm = TRUE), 1))
message("done")
