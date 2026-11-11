# Step 2b: cell-type annotation (marker-based, res 0.8 + neuronal sub-clustering) and verified region labels.
# Evidence for each label: results/annotation/markers_dotplot_res0.8.png, cluster_markers_res0.8.csv,
# and the neuronal sub-cluster marker table printed by this script.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
set.seed(42)
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
out <- file.path(proj, "results/annotation")
so <- readRDS(file.path(proj, "data/seurat_working_hvg.rds"))

# ---- region: file-name labels verified (GEO sample titles swap Hippocampus/Striatum) ----
# Striatum files: highest PPP1R1B+/PENK+ fraction in every animal + ependymal/SVZ (FOXJ1+, DLX+ neuroblasts);
# Hippocampus files: PROX1+/NEUROD2+ granule-like neurons, choroid plexus (TTR). See region_verification.csv.
so$region <- factor(so$region_file, levels = c("Amygdala", "Hippocampus", "Striatum"))
so$region_GEO_title <- unname(c(Amygdala = "Amygdala", Hippocampus = "Striatum", Striatum = "Hippocampus")[as.character(so$region)])
so$condition <- factor(so$condition, levels = c("Control", "ZIKV"))

fine_map <- c(
  "0" = "Microglia", "3" = "Microglia_activated (MHC-II+/CCL3/4+)", "19" = "Microglia_FOSB+",
  "29" = "Microglia_cycling", "34" = "Microglia_neuronal-transcript+",
  "26" = "Border-associated macrophages", "32" = "Monocytes", "31" = "T/NK cells", "36" = "B cells",
  "1" = "Astrocytes_protoplasmic", "8" = "Astrocytes_protoplasmic", "11" = "Astrocytes_SLC6A11+",
  "6" = "Astrocytes_low-UMI", "4" = "Astrocytes_GFAP-hi", "13" = "Astrocytes_immature (NES+/VIM+)",
  "35" = "Astrocytes_neuronal-transcript+",
  "2" = "OPC", "7" = "OPC", "21" = "OPC_cycling", "17" = "Cycling progenitors",
  "25" = "Neural progenitors (ASCL1+/DLL3+)",
  "9" = "COP/newly formed oligodendrocytes", "5" = "Oligodendrocytes", "14" = "Oligodendrocytes",
  "16" = "Oligodendrocytes", "24" = "Oligodendrocytes", "23" = "Oligodendrocytes_M5-hippocampus-specific",
  "12" = "Endothelial", "28" = "Pericytes", "15" = "Smooth muscle cells", "30" = "Fibroblasts",
  "27" = "Ependymal",
  "20" = "Low quality/doublet", "33" = "Low quality/doublet", "37" = "Low quality/doublet", "38" = "Low quality/doublet",
  "10" = "NEURON", "18" = "NEURON", "22" = "NEURON")
so$celltype_fine <- unname(fine_map[as.character(so$RNA_snn_res.0.8)])
stopifnot(!any(is.na(so$celltype_fine)))

# ---- neuronal sub-clustering ----
neu <- subset(so, celltype_fine == "NEURON")
neu <- FindNeighbors(neu, reduction = "harmony", dims = 1:20, verbose = FALSE) |> FindClusters(resolution = 0.5, verbose = FALSE)
neu_map <- c("0" = "Immature inhibitory neurons (DLX+/SP9+)", "2" = "Immature inhibitory neurons (DLX+/SP9+)",
             "7" = "Immature inhibitory neurons (DLX+/SP9+)",
             "1" = "Immature excitatory neurons (NEUROD2+/PROX1+)", "3" = "Immature excitatory neurons (NEUROD2+/PROX1+)",
             "4" = "Immature excitatory neurons (NEUROD2+/PROX1+)",
             "5" = "Neurons_PENK+", "6" = "Inhibitory interneurons", "8" = "Inhibitory interneurons",
             "9" = "Low quality/doublet", "10" = "Low quality/doublet")
g <- c("SNAP25","SYT1","SLC17A7","SLC17A6","NEUROD2","PROX1","GAD1","GAD2","SLC32A1","DLX2","SP9","SST","VIP","PENK","PPP1R1B","PDGFRA","MBP","AQP4","P2RY12")
print(round(t(as.matrix(AverageExpression(neu, features = g, layer = "data")$RNA)), 1))
ctf <- setNames(as.character(so$celltype_fine), colnames(so))
ctf[colnames(neu)] <- unname(neu_map[as.character(neu$seurat_clusters)])
stopifnot(!any(is.na(ctf)), !any(ctf == "NEURON"))
so$celltype_fine <- unname(ctf)

broad <- function(x) dplyr::case_when(
  grepl("^Microglia", x) ~ "Microglia",
  grepl("^Astrocytes", x) ~ "Astrocytes",
  x %in% c("OPC", "OPC_cycling") ~ "OPC",
  grepl("^Oligodendrocytes", x) ~ "Oligodendrocytes",
  TRUE ~ x)
so$celltype <- broad(unname(ctf))
ct_order <- c("Astrocytes", "Microglia", "Border-associated macrophages", "Monocytes", "T/NK cells", "B cells",
              "OPC", "Cycling progenitors", "COP/newly formed oligodendrocytes", "Oligodendrocytes",
              "Neural progenitors (ASCL1+/DLL3+)", "Immature inhibitory neurons (DLX+/SP9+)",
              "Immature excitatory neurons (NEUROD2+/PROX1+)", "Inhibitory interneurons", "Neurons_PENK+",
              "Endothelial", "Pericytes", "Smooth muscle cells", "Fibroblasts", "Ependymal", "Low quality/doublet")
stopifnot(all(unique(so$celltype) %in% ct_order))
so$celltype <- factor(so$celltype, levels = ct_order)
so$celltype_fine <- factor(so$celltype_fine)
Idents(so) <- "celltype"

# ---- outputs ----
write.csv(as.data.frame.matrix(table(so$celltype, so$sample)), file.path(out, "celltype_counts_per_sample.csv"))
write.csv(as.data.frame.matrix(table(so$celltype_fine, so$sample)), file.path(out, "celltype_fine_counts_per_sample.csv"))
write.csv(data.frame(cell = colnames(so), so@meta.data[, c("sample","animal","condition","region","region_GEO_title",
          "RNA_snn_res.0.8","celltype","celltype_fine")]), file.path(out, "cell_annotations.csv"), row.names = FALSE)
dp_genes <- list(Astro = c("AQP4","GFAP","SLC1A2","ALDH1L1"), Micro = c("P2RY12","CX3CR1","CD74"),
                 Mac = c("LYVE1","MRC1","CD163"), Mono = c("FCN1"), Lymph = c("CD3E","CD79B"),
                 OPC = c("PDGFRA","OLIG2"), Cyc = c("MKI67","TOP2A"), COP = c("BCAS1","GPR17"), Oligo = c("MOG","MBP"),
                 Prog = c("ASCL1","DLL3"), Neu = c("SNAP25","DLX2","SP9","NEUROD2","PROX1","GAD1","SST","VIP","PENK"),
                 Vasc = c("CLDN5","RGS5","ACTA2","COL1A1"), Epend = c("FOXJ1"), Target = c("DMD","DAG1"))
dp_genes <- lapply(dp_genes, intersect, rownames(so))
p <- DotPlot(subset(so, celltype != "Low quality/doublet"), features = dp_genes) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 8), strip.text = element_text(size = 8)) + labs(x = NULL, y = NULL)
ggsave(file.path(out, "annotation_marker_dotplot.png"), p, width = 18, height = 8, dpi = 200)
print(table(so$celltype, so$condition))
saveRDS(so, file.path(proj, "data/seurat_working_hvg_annotated.rds"))
message("done")
