# Step 2a: marker exploration for annotation + verification of region labels (GEO titles vs file names disagree)
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(tidyr) })
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
out <- file.path(proj, "results/annotation"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
so <- readRDS(file.path(proj, "data/seurat_working_hvg.rds"))
Idents(so) <- "RNA_snn_res.0.8"

markers <- list(
  Neuron = c("RBFOX3", "SNAP25", "SYT1"),
  Excit = c("SLC17A7", "SLC17A6", "SATB2", "NEUROD6", "TBR1"),
  Inhib = c("GAD1", "GAD2", "SLC32A1", "SST", "PVALB", "VIP", "LAMP5", "LHX6", "ADARB2"),
  MSN = c("PPP1R1B", "DRD1", "DRD2", "PENK", "TAC1", "FOXP1", "RARB"),
  DG_granule = c("PROX1", "C1QL2", "CALB1"),
  Amyg = c("TSHZ1", "FOXP2", "ZFHX3", "MEIS2", "LHX9"),
  Astro = c("AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1"),
  OPC = c("PDGFRA", "CSPG4", "OLIG2"),
  Oligo = c("MBP", "MOG", "PLP1", "MOBP", "ST18"),
  COP = c("BCAS1", "GPR17", "TCF7L2"),
  Micro = c("CX3CR1", "P2RY12", "CSF1R", "C1QB", "PTPRC"),
  Endo_Mural = c("CLDN5", "FLT1", "PECAM1", "PDGFRB", "RGS5", "ACTA2"),
  Ependymal_CP = c("FOXJ1", "TTR", "CFAP299"),
  Fibro = c("COL1A2", "DCN", "LUM"),
  T_cell = c("CD3E", "CD2"),
  Prolif = c("MKI67", "TOP2A"))
markers <- lapply(markers, intersect, rownames(so))
p <- DotPlot(so, features = markers, cluster.idents = TRUE) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 7), strip.text = element_text(size = 7))
ggsave(file.path(out, "markers_dotplot_res0.8.png"), p, width = 26, height = 13, dpi = 150)

top <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, max.cells.per.ident = 500)
top10 <- top |> group_by(cluster) |> slice_max(avg_log2FC, n = 12)
write.csv(top, file.path(out, "cluster_markers_res0.8.csv"), row.names = FALSE)

# region verification: fraction of each sample's cells in MSN-like and DG-granule-like states
ms <- AddModuleScore(so, features = list(markers$MSN[1:4], c("PROX1", "C1QL2")), name = "sig")
v <- ms@meta.data |> group_by(sample, region_file) |>
  summarise(n = n(), MSN_frac = mean(sig1 > 0.5), DG_frac = mean(sig2 > 0.5),
            PPP1R1B = mean(FetchData(ms, "PPP1R1B")[cur_group_rows(), 1] > 0),
            PROX1 = mean(FetchData(ms, "PROX1")[cur_group_rows(), 1] > 0), .groups = "drop")
print(as.data.frame(v))
write.csv(v, file.path(out, "region_verification.csv"), row.names = FALSE)

cl <- table(so$RNA_snn_res.0.8, so$sample)
write.csv(as.data.frame.matrix(cl), file.path(out, "cluster_by_sample_res0.8.csv"))
qc <- so@meta.data |> group_by(RNA_snn_res.0.8) |>
  summarise(n = n(), umi = median(nCount_RNA), genes = median(nFeature_RNA), mt = median(percent.mt),
            dbl = median(scDblFinder.score), n_samples = n_distinct(sample))
print(as.data.frame(qc))
for (cl_id in levels(so$RNA_snn_res.0.8)) {
  cat(cl_id, ":", paste(head(top10$gene[top10$cluster == cl_id], 12), collapse = ", "), "\n")
}
p <- DimPlot(so, label = TRUE) + NoLegend()
ggsave(file.path(out, "umap_res0.8.png"), p, width = 8, height = 7, dpi = 150)
