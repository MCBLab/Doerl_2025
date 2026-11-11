# Step 5: data for the report linking these results to Doerl et al. 2025 (bioRxiv 10.1101/2025.10.26.684625,
# "Maturation of the glymphatic system confers innate resistance of the brain to Zika virus infection").
# Extends the AQP4/DMD/DAG1 analysis to the rest of the astrocyte endfoot / dystrophin-glycoprotein complex (DGC),
# its basement-membrane ligands, vascular tone (incl. AVPR1A = V1aR, the receptor targeted in the preprint),
# BBB genes and interferon-stimulated genes. Same pseudobulk edgeR model as 03; no cell-level statistics.
# Light: uses only the pseudobulk counts and the genome-wide tables written by 03 (Seurat object not loaded).
suppressPackageStartupMessages({ library(edgeR); library(dplyr); library(jsonlite) })
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
tab_dir <- file.path(proj, "results/tables")
MIN_CELLS <- 10
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

panel <- list(
  "Endfoot / DGC (astrocyte side)" = c("AQP4", "KCNJ10", "DMD", "DAG1", "SNTA1", "SNTB1", "DTNA", "DTNB", "UTRN",
                                       "MLC1", "HEPACAM", "LARGE1", "POMT1"),
  "Basement-membrane ligands" = c("AGRN", "LAMA2", "LAMA4", "LAMB1", "LAMB2", "LAMC1", "HSPG2", "NID1", "COL4A1", "COL4A2"),
  "Vascular tone and BBB" = c("AVPR1A", "ACTA2", "MYH11", "PDGFRB", "CLDN5", "OCLN", "TJP1"),
  "Astrocyte reactivity" = c("GFAP", "VIM", "CD44"),
  "Interferon response" = c("ISG15", "IFI6", "IFI27", "MX1", "IFIT1", "IFIT3", "STAT1", "OAS1"))
panel_df <- data.frame(group = rep(names(panel), lengths(panel)), gene = unlist(panel), row.names = NULL)
pb <- readRDS(file.path(proj, "data/pseudobulk_animal_region_celltype.rds"))
stopifnot(all(panel_df$gene %in% rownames(pb$counts_broad)))

# ---- where each panel gene is expressed: mean log2 CPM in Control pseudobulks, by cell type ----
lib <- colSums(pb$counts_broad)
cpm_all <- t(t(pb$counts_broad[unique(panel_df$gene), ]) / lib * 1e6)
expr <- bind_rows(lapply(split(seq_len(ncol(cpm_all)), pb$info_broad$celltype), function(j) {
  j <- j[pb$info_broad$condition[j] == "Control" & pb$info_broad$n_cells[j] >= MIN_CELLS]
  if (length(j) < 3) return(NULL)
  data.frame(celltype = as.character(pb$info_broad$celltype[j[1]]), gene = rownames(cpm_all),
             log2CPM = log2(rowMeans(cpm_all[, j, drop = FALSE]) + 1), n_pb = length(j))
}))

# ---- ambient-RNA / doublet contamination check (scripts/contamination.R, shared with 03b) ----
source(file.path(proj, "scripts/contamination.R"))
P_broad <- pooled_cpm(pb$counts_broad, with(pb$info_broad, paste(region, condition, celltype, sep = "|")), panel_df$gene)
P_fine <- pooled_cpm(pb$counts_fine, with(pb$info_fine, paste(region, condition, celltype_fine, sep = "|")), panel_df$gene)
flag_contam <- function(d, P, ct, ct_broad, rg) {
  d <- left_join(d, contamination(panel_df$gene, P, ct, ct_broad, rg, P_broad), by = "gene")
  d$ambient_prone <- d$tested & d$ambient_prone
  d
}

# ---- broad cell types: panel genes from the genome-wide edgeR tables written by 03 ----
tested <- pb$info_broad |> filter(n_cells >= MIN_CELLS) |> count(region, celltype, condition) |>
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) |> filter(Control >= 2, ZIKV >= 2)
broad <- bind_rows(lapply(seq_len(nrow(tested)), function(k) {
  rg <- as.character(tested$region[k]); ct <- as.character(tested$celltype[k])
  tt <- read.csv(file.path(tab_dir, "genome_wide", sprintf("edgeR_pseudobulk_%s_%s.csv", rg, safe(ct))))
  d <- left_join(panel_df, tt[, c("gene", "logFC", "logCPM", "PValue", "FDR")], by = "gene") |>
    mutate(region = rg, celltype = ct, tested = !is.na(PValue))
  flag_contam(d, P_broad, ct, ct, rg)
}))
write.csv(broad, file.path(tab_dir, "glymphatic_panel_pseudobulk_edgeR_broad.csv"), row.names = FALSE)

# ---- astrocyte subtypes (fine level): rerun the same model, keep panel genes ----
run_fine <- function(rg, ct) {
  info <- pb$info_fine |> mutate(i = row_number()) |>
    filter(region == rg, celltype_fine == ct, n_cells >= MIN_CELLS)
  if (sum(info$condition == "Control") < 2 || sum(info$condition == "ZIKV") < 2) return(NULL)
  y <- DGEList(pb$counts_fine[, info$i], group = factor(info$condition, levels = c("Control", "ZIKV")))
  y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
  design <- model.matrix(~ group, data = y$samples)
  fit <- glmQLFit(estimateDisp(y, design), design, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table; tt$gene <- rownames(tt)
  d <- left_join(panel_df, tt[, c("gene", "logFC", "logCPM", "PValue", "FDR")], by = "gene") |>
    mutate(region = rg, celltype_fine = ct, tested = !is.na(PValue))
  flag_contam(d, P_fine, ct, "Astrocytes", rg)
}
astro_sub <- c("Astrocytes_GFAP-hi", "Astrocytes_protoplasmic")
fine <- bind_rows(lapply(levels(pb$info_fine$region), function(rg) bind_rows(lapply(astro_sub, run_fine, rg = rg))))
write.csv(fine, file.path(tab_dir, "glymphatic_panel_pseudobulk_edgeR_astrocyte_subtypes.csv"), row.names = FALSE)

# ---- per-animal CPM for the key endfoot genes in astrocytes and GFAP-hi astrocytes ----
key <- c("AQP4", "DMD", "DAG1", "KCNJ10", "UTRN", "NID1", "LAMB2", "HSPG2")
pa <- function(mat, info, lvl) {
  ok <- info$n_cells >= MIN_CELLS
  y <- normLibSizes(DGEList(mat[, ok])); c <- cpm(y)[key, , drop = FALSE]
  data.frame(info[ok, c("animal", "condition", "region", "n_cells")], level = as.character(info[[lvl]][ok]),
             t(c), check.names = FALSE, row.names = NULL)
}
cpm_astro <- bind_rows(
  pa(pb$counts_broad[, pb$info_broad$celltype == "Astrocytes"], pb$info_broad[pb$info_broad$celltype == "Astrocytes", ], "celltype"),
  pa(pb$counts_fine[, pb$info_fine$celltype_fine == "Astrocytes_GFAP-hi"],
     pb$info_fine[pb$info_fine$celltype_fine == "Astrocytes_GFAP-hi", ], "celltype_fine"))
write.csv(cpm_astro, file.path(tab_dir, "glymphatic_key_genes_CPM_per_animal_astrocytes.csv"), row.names = FALSE)

# ---- composition: cell-type proportions per sample (QC'd singlets, doublet clusters removed) ----
comp <- pb$info_broad |> group_by(animal, condition, region, celltype) |> summarise(n = sum(n_cells), .groups = "drop") |>
  group_by(animal, condition, region) |> mutate(total = sum(n), pct = 100 * n / total) |> ungroup()
write.csv(comp, file.path(tab_dir, "celltype_proportions_per_sample.csv"), row.names = FALSE)

# ---- export for the HTML report ----
num <- function(d) mutate(d, across(where(is.numeric), ~ signif(.x, 4)))
out <- list(panel = panel, broad = num(broad), fine = num(fine), expr = num(expr),
            cpm_astro = num(cpm_astro), comp = num(comp))
write_json(out, file.path(proj, "data/report_data.json"), dataframe = "rows", na = "null", digits = NA)

# console summary
show <- function(x) base::print(tibble::as_tibble(x), n = 200, width = 200)
show(broad |> filter(tested, PValue < 0.05, celltype %in% c("Astrocytes", "Endothelial", "Pericytes", "Smooth muscle cells")) |>
       select(region, celltype, gene, logFC, PValue, FDR, log2CPM_ctrl, log2FC_pooled, log2FC_corrected, ambient_prone) |> arrange(celltype, region, PValue))
show(fine |> filter(tested, PValue < 0.05) |> select(region, celltype_fine, gene, logFC, PValue, FDR, log2FC_pooled, log2FC_corrected, ambient_prone) |>
       arrange(celltype_fine, region, PValue))
show(broad |> filter(gene == "AVPR1A", tested) |> select(region, celltype, logFC, PValue, FDR))
message("done")
