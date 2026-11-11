# Step 1b: integrate the 15 per-sample QC'd objects (Harmony), cluster, UMAP.
# Memory-lean: the full merged matrix (~179k cells, ~500M non-zeros) needs ~6 GB per copy, so integration is
# run on a working object holding only the integration HVGs + marker/target genes. Log-normalisation is done
# per sample on the FULL transcriptome before subsetting, so normalised values are identical to a full object.
suppressPackageStartupMessages({ library(Seurat); library(harmony); library(dplyr); library(Matrix) })
set.seed(42)
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
ckpt_dir <- file.path(proj, "data/per_sample_ckpt")
out_dir <- file.path(proj, "results/qc"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- sort(list.files(ckpt_dir, full.names = TRUE))

targets <- c("AQP4", "DMD", "DAG1")
marker_panel <- c("RBFOX3","SNAP25","SYT1","SLC17A7","SLC17A6","SATB2","NEUROD6","TBR1","GAD1","GAD2","SLC32A1",
  "SST","PVALB","VIP","LAMP5","LHX6","ADARB2","PPP1R1B","DRD1","DRD2","PENK","TAC1","FOXP1","RARB","PROX1","C1QL2",
  "CALB1","TSHZ1","FOXP2","ZFHX3","MEIS2","LHX9","GFAP","SLC1A2","SLC1A3","ALDH1L1","GJA1","PDGFRA","CSPG4","OLIG2",
  "MBP","MOG","PLP1","MOBP","ST18","BCAS1","GPR17","TCF7L2","CX3CR1","P2RY12","CSF1R","C1QB","PTPRC","CLDN5","FLT1",
  "PECAM1","PDGFRB","RGS5","ACTA2","FOXJ1","TTR","CFAP299","COL1A2","DCN","LUM","CD3E","CD2","MKI67","TOP2A",
  "CUX2","RORB","FEZF2","THEMIS","SLA","NR4A2","CCK","CNR1","RELN","NPY","CHAT","TH","SNTA1","UTRN","DTNA","KCNJ10")

pca_ckpt <- file.path(proj, "data/seurat_working_hvg_pca.rds")
if (!file.exists(pca_ckpt)) {
# ---- pass 1: per-sample HVGs, QC table, metadata ----
hvg_rank <- list(); qc_tab <- list(); mds <- list()
for (f in files) {
  r <- readRDS(f); s <- r$so; nm <- s$sample[1]
  s <- FindVariableFeatures(s, nfeatures = 3000, verbose = FALSE)
  hvg_rank[[nm]] <- setNames(seq_along(VariableFeatures(s)), VariableFeatures(s))
  qc_tab[[nm]] <- r$qc; mds[[nm]] <- s@meta.data
  if (is.null(genes_all <- get0("genes_all"))) genes_all <- rownames(s) else stopifnot(identical(genes_all, rownames(s)))
  rm(r, s); gc(verbose = FALSE)
}
qc_tab <- bind_rows(qc_tab); print(qc_tab)
write.csv(qc_tab, file.path(out_dir, "qc_summary_per_sample.csv"), row.names = FALSE)
# SelectIntegrationFeatures-like ranking: # samples where variable, then median rank
all_hvg <- unique(unlist(lapply(hvg_rank, names)))
n_var <- sapply(all_hvg, function(g) sum(sapply(hvg_rank, function(h) g %in% names(h))))
med_rank <- sapply(all_hvg, function(g) median(unlist(lapply(hvg_rank, function(h) h[g])), na.rm = TRUE))
hvg <- all_hvg[order(-n_var, med_rank)][1:3000]
keep_genes <- union(hvg, intersect(c(targets, marker_panel), genes_all))
message(length(keep_genes), " genes kept in working object")

# ---- pass 2: normalise each sample on full transcriptome, keep subset ----
cnt_l <- list(); dat_l <- list()
for (f in files) {
  s <- readRDS(f)$so
  s <- NormalizeData(s, verbose = FALSE)
  cnt_l[[s$sample[1]]] <- LayerData(s, layer = "counts")[keep_genes, ]
  dat_l[[s$sample[1]]] <- LayerData(s, layer = "data")[keep_genes, ]
  rm(s); gc(verbose = FALSE)
}
cnt <- do.call(cbind, cnt_l); dat <- do.call(cbind, dat_l); rm(cnt_l, dat_l)
md <- bind_rows(mds); md <- md[colnames(cnt), ]
md$nCount_RNA_full <- md$nCount_RNA; md$nFeature_RNA_full <- md$nFeature_RNA   # full-transcriptome QC metrics
so <- CreateSeuratObject(counts = cnt, meta.data = md[, setdiff(colnames(md), c("orig.ident","nCount_RNA","nFeature_RNA"))])
so <- SetAssayData(so, layer = "data", new.data = dat); rm(cnt, dat); gc()
so$nCount_RNA <- so$nCount_RNA_full; so$nFeature_RNA <- so$nFeature_RNA_full
so$nCount_RNA_full <- NULL; so$nFeature_RNA_full <- NULL
VariableFeatures(so) <- hvg
# PCA on the sparse log-normalised HVG matrix with implicit centring/scaling (irlba) — avoids the dense
# 3000 x 179k scale.data matrix (~4.3 GB, plus copies inside RunPCA). Equivalent to ScaleData+RunPCA
# except that z-scores are not clipped at 10.
Xt <- Matrix::t(LayerData(so, layer = "data")[hvg, ])
mu <- Matrix::colMeans(Xt)
sdv <- sqrt((Matrix::colSums(Xt^2) - nrow(Xt) * mu^2) / (nrow(Xt) - 1))
pc <- irlba::irlba(Xt, nv = 50, center = mu, scale = sdv)
rm(Xt); gc()
emb <- pc$u %*% diag(pc$d); dimnames(emb) <- list(colnames(so), paste0("PC_", 1:50))
lds <- pc$v; dimnames(lds) <- list(hvg, paste0("PC_", 1:50))
so[["pca"]] <- CreateDimReducObject(embeddings = emb, loadings = lds, key = "PC_", assay = "RNA",
                                    stdev = pc$d / sqrt(ncol(so) - 1))
rm(pc, emb, lds); gc()
saveRDS(so, pca_ckpt)   # checkpoint
} else so <- readRDS(pca_ckpt)
so <- RunHarmony(so, group.by.vars = "sample", dims.use = 1:30, project.dim = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = 1:30)
so <- FindNeighbors(so, reduction = "harmony", dims = 1:30) |> FindClusters(resolution = c(0.3, 0.8))
so@misc$note <- "Working object: integration HVGs + marker/target genes only; data layer log-normalised on full transcriptome per sample."
saveRDS(so, file.path(proj, "data/seurat_working_hvg.rds"))
message("done")
