# Step 8: violin plots of the endfoot-anchor genes in astrocytes, ZIKV vs Control.
# Genes: DAG1, DMD, UTRN, LARGE1, SNTA1, DTNA -- the dystrophin-glycoprotein complex that holds AQP4
# at the perivascular endfoot.
#
# These transcripts are sparsely detected per cell, so a violin of all cells is dominated by the zero
# mass. The figure therefore shows the full per-cell distribution (honest about the zeros), the percent
# of cells with any detection, and -- as the points -- the per-animal means, which are the actual unit
# of replication. The statistic annotated on each panel is the MAST hurdle test with the animal as a
# random effect (script 07), not a cell-level P-value.
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(dplyr) })
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
fig_dir <- file.path(proj, "results/figures"); tab_dir <- file.path(proj, "results/tables")

ANCHOR <- c("DAG1", "DMD", "UTRN", "LARGE1", "SNTA1", "DTNA")
CTRL <- "#2a78d6"; ZIKV <- "#eb6834"
pe <- readRDS(file.path(proj, "data/panel_cell_expr.rds"))
counts <- pe$counts; md <- pe$md; rm(pe)
md$condition <- factor(md$condition, levels = c("Control", "ZIKV"))
lognorm <- log2(t(t(counts[ANCHOR, , drop = FALSE]) / md$nCount_RNA) * 1e4 + 1)

theme_set(theme_minimal(base_size = 11))
th <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(colour = "grey91", linewidth = .35),
  strip.text = element_text(face = "bold", size = 11, hjust = 0),
  plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(colour = "grey35", size = 9.5, lineheight = 1.2),
  plot.caption = element_text(colour = "grey40", size = 8, hjust = 0, lineheight = 1.2),
  legend.position = "top", legend.justification = "left",
  axis.title.x = element_blank())

make_panel <- function(cells, label, stats_file, stat_key) {
  long <- bind_rows(lapply(ANCHOR, function(g) data.frame(
    gene = g, expr = lognorm[g, cells], detected = counts[g, cells] > 0,
    condition = md$condition[cells], animal = md$animal[cells], region = md$region[cells])))
  long$gene <- factor(long$gene, levels = ANCHOR)

  # every layer needs the same factor levels, or the facets fall back to alphabetical order
  fac <- function(d) { d$gene <- factor(d$gene, levels = ANCHOR); d }
  per_animal <- long |> group_by(gene, condition, animal, region) |>
    summarise(m = mean(expr), .groups = "drop") |> fac()
  pct <- long |> group_by(gene, condition) |>
    summarise(pct = 100 * mean(detected), n = n(), .groups = "drop")
  ymax <- long |> group_by(gene) |> summarise(top = quantile(expr, .999) * 1.3 + .15, .groups = "drop")
  pct <- left_join(pct, ymax, by = "gene") |> fac()

  st <- read.csv(file.path(tab_dir, stats_file)) |> filter(.data[[stat_key]] == label, gene %in% ANCHOR)
  lab <- ymax |> left_join(st, by = "gene") |> fac() |>
    mutate(txt = ifelse(is.na(PValue), "not testable",
                 sprintf("log2FC %+.2f · P %s%s", logFC_MAST,
                         ifelse(PValue < 0.001, format(PValue, digits = 2, scientific = TRUE), sprintf("%.3f", PValue)),
                         ifelse(!is.na(FDR) & FDR < 0.05, " · FDR<0.05", ""))))

  ggplot(long, aes(condition, expr, fill = condition)) +
    geom_violin(scale = "width", width = .78, linewidth = .3, colour = "grey45", trim = TRUE,
                alpha = .55, show.legend = TRUE) +
    geom_point(data = per_animal, aes(x = condition, y = m, colour = condition, shape = region),
               position = position_jitter(width = .13, height = 0, seed = 1),
               size = 2.1, stroke = .9, fill = "white", inherit.aes = FALSE) +
    geom_text(data = pct, aes(x = condition, y = top * .04, label = sprintf("%.0f%%", pct)),
              inherit.aes = FALSE, size = 2.9, colour = "grey30", vjust = 1.4) +
    geom_text(data = lab, aes(x = 1.5, y = top, label = txt), inherit.aes = FALSE,
              size = 2.9, colour = "grey25") +
    facet_wrap(~ gene, nrow = 2, scales = "free_y") +
    scale_fill_manual(values = c(Control = CTRL, ZIKV = ZIKV), name = NULL) +
    scale_colour_manual(values = c(Control = CTRL, ZIKV = ZIKV), guide = "none") +
    scale_shape_manual(values = c(Amygdala = 21, Hippocampus = 22, Striatum = 24), name = "Region (per-animal mean)") +
    scale_y_continuous(expand = expansion(mult = c(.04, .14))) +
    labs(y = expression(log[2]*"(counts per 10,000 + 1)")) + th
}

cap <- paste(
  "Each violin is the full per-cell distribution, width-scaled; percentages under each violin are the share of cells",
  "with any detection. Points are per-animal means (3 control, 2 ZIKV animals × 3 regions), which are the unit of",
  "replication. The annotated statistic is the MAST hurdle test with the animal as a random effect and region as a",
  "covariate (script 07); cell-level P-values are not used. n = 3 vs 2 animals — hypothesis-generating.", sep = "\n")

p1 <- make_panel(which(md$celltype == "Astrocytes"), "Astrocytes",
                 "mast_endfoot_anchor_animalRE_by_celltype.csv", "celltype") +
  labs(title = "Endfoot anchor complex in astrocytes, ZIKV vs Control",
       subtitle = paste0("Infant rhesus macaque, 14 dpi · ",
                         format(sum(md$celltype == "Astrocytes"), big.mark = ","),
                         " astrocytes pooled across amygdala, hippocampus and striatum"),
       caption = cap)
ggsave(file.path(fig_dir, "08_violin_endfoot_astrocytes.png"), p1, width = 11, height = 7.2, dpi = 200, bg = "white")

p2 <- make_panel(which(md$celltype_fine == "Astrocytes_GFAP-hi"), "Astrocytes_GFAP-hi",
                 "mast_endfoot_anchor_animalRE_astrocyte_subtypes.csv", "celltype_fine") +
  labs(title = "Endfoot anchor complex in GFAP-high astrocytes, ZIKV vs Control",
       subtitle = paste0("The subtype where the pseudobulk effect is strongest · ",
                         format(sum(md$celltype_fine == "Astrocytes_GFAP-hi"), big.mark = ","), " cells"),
       caption = cap)
ggsave(file.path(fig_dir, "09_violin_endfoot_astrocytes_GFAPhi.png"), p2, width = 11, height = 7.2, dpi = 200, bg = "white")

message("wrote 08_violin_endfoot_astrocytes.png and 09_violin_endfoot_astrocytes_GFAPhi.png")
