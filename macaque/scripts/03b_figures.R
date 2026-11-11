# Step 3b: figures + ambient-RNA diagnostic for the AQP4 / DMD / DAG1 ZIKV-vs-Control comparison.
# Reads results produced by 03_de_target_genes.R (no re-computation of DE).
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
fig_dir <- file.path(proj, "results/figures"); tab_dir <- file.path(proj, "results/tables")
so <- readRDS(file.path(proj, "data/seurat_working_hvg_annotated.rds"))
genes <- c("AQP4", "DMD", "DAG1")
cond_cols <- c(Control = "#4C78A8", ZIKV = "#E45756")
MIN_CELLS <- 10
res <- so@misc$target_gene_DE$pseudobulk_edgeR
cpm_df <- so@misc$target_gene_DE$pseudobulk_cpm
ct_levels <- levels(so$celltype)
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                          strip.background = element_rect(fill = "grey95")))

# ---------------- ambient-RNA diagnostic ----------------
# AQP4 (and partly DMD) are astrocyte-enriched. If ZIKV libraries carry more astrocyte-derived ambient RNA,
# non-astrocyte pseudobulks will show spurious "increases" of these genes together with other astrocyte-only
# genes. Index = mean ZIKV-vs-Control logFC of 5 astrocyte-specific genes in the same region x cell type.
astro_markers <- c("GFAP", "SLC1A2", "SLC1A3", "GJA1", "ALDOC")
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
res$astro_ambient_logFC <- NA_real_
for (k in which(!is.na(res$PValue) & res$celltype != "Astrocytes")) {
  f <- file.path(tab_dir, "genome_wide", sprintf("edgeR_pseudobulk_%s_%s.csv", res$region[k], safe(res$celltype[k])))
  tt <- read.csv(f)
  res$astro_ambient_logFC[k] <- mean(tt$logFC[tt$gene %in% astro_markers])
}
# Flag (†): the quantitative contamination check of scripts/contamination.R (all source cell types, not only
# astrocytes; shared with 05). The astrocyte index above is kept for the diagnostic figure.
source(file.path(proj, "scripts/contamination.R"))
pb <- readRDS(file.path(proj, "data/pseudobulk_animal_region_celltype.rds"))
P_broad <- pooled_cpm(pb$counts_broad, with(pb$info_broad, paste(region, condition, celltype, sep = "|")), genes)
cc <- bind_rows(lapply(split(res, list(res$region, res$celltype), drop = TRUE), function(d) {
  rg <- as.character(d$region[1]); ct <- as.character(d$celltype[1])
  if (!all(paste(rg, c("Control", "ZIKV"), ct, sep = "|") %in% colnames(P_broad))) return(NULL)
  data.frame(region = rg, celltype = ct, contamination(genes, P_broad, ct, ct, rg, P_broad))
}))
rm(pb, P_broad)
res <- res |> select(-any_of(c("ambient_caution", setdiff(colnames(cc), c("region", "celltype", "gene"))))) |>
  left_join(cc, by = c("region", "celltype", "gene")) |> mutate(ambient_prone = !is.na(PValue) & ambient_prone)
write.csv(res, file.path(tab_dir, "target_genes_pseudobulk_edgeR_ZIKV_vs_Control.csv"), row.names = FALSE)
so@misc$target_gene_DE$pseudobulk_edgeR <- res

amb <- res |> filter(!is.na(logFC), !is.na(astro_ambient_logFC))
r_lab <- amb |> group_by(gene) |> summarise(r = cor(logFC, astro_ambient_logFC), .groups = "drop") |>
  mutate(lab = sprintf("r = %.2f", r))
p <- ggplot(amb, aes(astro_ambient_logFC, logFC)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = region), size = 2.2, alpha = 0.85) +
  geom_text(data = r_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.2, vjust = 1.5, inherit.aes = FALSE) +
  facet_wrap(~ gene) +
  labs(x = "Astrocyte-specific gene logFC (mean of GFAP, SLC1A2, SLC1A3, GJA1, ALDOC)",
       y = "Target gene logFC, ZIKV vs Control",
       title = "Ambient-RNA diagnostic in NON-astrocyte cell types",
       subtitle = "AQP4/DMD shifts that track astrocyte-only genes reflect astrocyte-derived ambient RNA, not regulation in that cell type")
ggsave(file.path(fig_dir, "07_ambient_RNA_diagnostic.png"), p, width = 12, height = 4.8, dpi = 200)

# ---------------- 1. annotated UMAPs ----------------
p <- DimPlot(so, group.by = "celltype", label = TRUE, repel = TRUE, label.size = 3) + NoLegend() + ggtitle("Cell types") |
     DimPlot(so, group.by = "region", shuffle = TRUE) + ggtitle("Region") |
     DimPlot(so, group.by = "condition", shuffle = TRUE, cols = cond_cols) + ggtitle("Condition")
ggsave(file.path(fig_dir, "01_UMAP_celltype_region_condition.png"), p, width = 20, height = 6.5, dpi = 200)

# ---------------- 2. feature plots split by condition ----------------
p <- FeaturePlot(so, genes, split.by = "condition", order = TRUE, combine = TRUE, keep.scale = "feature") &
     theme(legend.position = "right")
ggsave(file.path(fig_dir, "02_UMAP_target_genes_by_condition.png"), p, width = 10, height = 12, dpi = 200)

# ---------------- 3. dot plot: celltype x condition, faceted by region ----------------
ed <- FetchData(so, vars = c(genes, "celltype", "condition", "region"), layer = "data")
avg <- ed |> pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(region, celltype, condition, gene) |>
  summarise(avg = mean(expm1(expr)), pct = 100 * mean(expr > 0), n = n(), .groups = "drop") |> filter(n >= 20)
avg <- avg |> group_by(gene) |> mutate(scaled = (log1p(avg) - mean(log1p(avg))) / sd(log1p(avg))) |> ungroup()
avg$celltype <- factor(avg$celltype, levels = rev(ct_levels)); avg$gene <- factor(avg$gene, levels = genes)
p <- ggplot(avg, aes(condition, celltype)) +
  geom_point(aes(size = pct, colour = scaled)) +
  facet_grid(gene ~ region) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey90", high = "#B40426", name = "Scaled\nlog mean") +
  scale_size_area(max_size = 6, name = "% cells\nexpressing") +
  labs(x = NULL, y = NULL, title = "AQP4 / DMD / DAG1 by cell type, region and condition (cell level)",
       subtitle = "Groups with >= 20 cells shown") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "03_dotplot_region_celltype_condition.png"), p, width = 10, height = 16, dpi = 200)

# ---------------- 4. pseudobulk CPM per animal (unit of replication), tested groups only ----------------
tested <- res |> filter(!is.na(PValue)) |> distinct(region, celltype)
cpm_long <- cpm_df |> pivot_longer(all_of(genes), names_to = "gene", values_to = "CPM") |>
  filter(n_cells >= MIN_CELLS) |> semi_join(tested, by = c("region", "celltype"))
cpm_long$celltype <- factor(cpm_long$celltype, levels = ct_levels)
lab <- res |> filter(!is.na(PValue)) |>
  mutate(label = paste0(case_when(FDR < 0.05 ~ "**", PValue < 0.05 ~ "*", TRUE ~ ""), ifelse(ambient_prone, "†", "")))
for (g in genes) {
  d <- filter(cpm_long, gene == g)
  l <- filter(lab, gene == g) |>
    inner_join(d |> group_by(region, celltype) |> summarise(y = max(log2(CPM + 1)), .groups = "drop"), by = c("region", "celltype"))
  l$celltype <- factor(l$celltype, levels = ct_levels)
  p <- ggplot(d, aes(celltype, log2(CPM + 1), fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.35, position = position_dodge(0.8), width = 0.7) +
    geom_point(aes(colour = condition), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.8) +
    geom_text(data = l, aes(celltype, y + 0.7, label = label), inherit.aes = FALSE, size = 5) +
    facet_wrap(~ region, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = cond_cols) + scale_colour_manual(values = cond_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = "log2(pseudobulk CPM + 1)", title = sprintf("%s: pseudobulk expression per animal", g),
         subtitle = paste("Each dot = one animal (Control n = 3, ZIKV n = 2); groups with >= 2 animals/arm and >= 10 cells.",
                          "edgeR QL: * nominal P < 0.05, ** genome-wide FDR < 0.05;",
                          "† = change not robust to removing the estimated ambient-RNA/doublet signal", sep = "\n")) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  ggsave(file.path(fig_dir, sprintf("04_pseudobulk_CPM_%s.png", g)), p, width = 11, height = 12, dpi = 200)
}

# ---------------- 5. log2FC heatmap (ZIKV vs Control) ----------------
hm <- res |> mutate(celltype = factor(celltype, levels = rev(ct_levels)), gene = factor(gene, levels = genes),
                    sig = paste0(case_when(FDR < 0.05 ~ "**", PValue < 0.05 ~ "*", TRUE ~ ""), ifelse(ambient_prone, "†", "")))
lim <- max(abs(hm$logFC), na.rm = TRUE)
p <- ggplot(hm, aes(gene, celltype, fill = logFC)) + geom_tile(colour = "white") +
  geom_text(aes(label = sig), size = 4.5, vjust = 0.75) +
  facet_wrap(~ region, nrow = 1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", limits = c(-lim, lim), na.value = "grey80",
                       name = "log2FC\nZIKV vs Ctrl") +
  labs(x = NULL, y = NULL, title = "AQP4 / DMD / DAG1: pseudobulk log2 fold-change, ZIKV vs Control",
       caption = paste("edgeR quasi-likelihood on animal-level pseudobulk (Control n = 3, ZIKV n = 2).",
                       "* nominal P < 0.05; ** genome-wide FDR < 0.05; † change not robust to removing the estimated ambient-RNA/doublet signal.",
                       "Grey = gene too lowly expressed to test; blank row = < 2 animals per arm with >= 10 cells.", sep = "\n")) +
  theme(panel.grid = element_blank(), plot.caption = element_text(hjust = 0, size = 9))
ggsave(file.path(fig_dir, "05_heatmap_log2FC_ZIKV_vs_Control.png"), p, width = 11, height = 7.5, dpi = 200)

# ---------------- 6. violin per region (cell level, descriptive) ----------------
for (rg in levels(so$region)) {
  s <- subset(so, region == rg)
  keep_ct <- names(which(apply(table(s$celltype, s$condition), 1, min) >= 20))
  s <- subset(s, celltype %in% keep_ct); s$celltype <- droplevels(s$celltype)
  p <- VlnPlot(s, genes, group.by = "celltype", split.by = "condition", cols = cond_cols, pt.size = 0,
               ncol = 1, alpha = 0.1) & theme(axis.title.x = element_blank())
  p <- p + plot_annotation(title = sprintf("%s: cell-level expression (descriptive; statistics in pseudobulk figures)", rg))
  ggsave(file.path(fig_dir, sprintf("06_violin_%s.png", rg)), p, width = 13, height = 12, dpi = 200)
}

saveRDS(so, file.path(proj, "data/seurat_working_hvg_annotated.rds"))
message("done")
