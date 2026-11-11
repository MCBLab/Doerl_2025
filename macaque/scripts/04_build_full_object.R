# Step 4: assemble the deliverable — full-transcriptome annotated Seurat object.
# Memory-lean: counts are streamed from per-sample checkpoints into pre-allocated slot vectors; the normalised
# "data" layer shares i/p with counts (only its x costs memory). The Seurat object is created around empty
# placeholder layers and the real matrices are put in the layer slots directly — CreateSeuratObject() on the
# real matrix makes ~2.7x transient copies (>16 GB here). The object is validated before saving.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
wk <- readRDS(file.path(proj, "data/seurat_working_hvg_annotated.rds"))   # low-quality/doublet clusters already removed
# keep only the light parts of the working object; its expression layers (~2.4 GB) are not needed
keep <- list(meta = wk@meta.data, red = lapply(setNames(Reductions(wk), Reductions(wk)), function(r) wk[[r]]),
             graphs = lapply(setNames(Graphs(wk), Graphs(wk)), function(g) wk[[g]]),
             hvg = VariableFeatures(wk), misc = wk@misc)
cells <- colnames(wk); rm(wk); gc()
files <- sort(list.files(file.path(proj, "data/per_sample_ckpt"), full.names = TRUE))

# pass 1: sizes
nnz <- 0; genes <- NULL; ord <- character(0)
for (f in files) {
  m <- LayerData(readRDS(f)$so, layer = "counts"); m <- m[, colnames(m) %in% cells, drop = FALSE]
  nnz <- nnz + length(m@x); ord <- c(ord, colnames(m))
  if (is.null(genes)) genes <- rownames(m) else stopifnot(identical(genes, rownames(m)))
  rm(m); gc(verbose = FALSE)
}
stopifnot(setequal(ord, cells))
message("cells: ", length(ord), "  non-zeros: ", nnz)
x <- numeric(nnz); dx <- numeric(nnz); i <- integer(nnz); p <- integer(length(ord) + 1)
# pass 2: fill
pos <- 0; cpos <- 0L
for (f in files) {
  m <- LayerData(readRDS(f)$so, layer = "counts"); m <- m[, colnames(m) %in% cells, drop = FALSE]
  k <- length(m@x); idx <- pos + seq_len(k)
  x[idx] <- m@x; i[idx] <- m@i
  lib <- colSums(m)
  dx[idx] <- log1p(m@x / rep(lib, diff(m@p)) * 1e4)        # Seurat LogNormalize, scale.factor 1e4
  p[cpos + seq_len(ncol(m)) + 1L] <- as.integer(m@p[-1] + pos)
  pos <- pos + k; cpos <- cpos + ncol(m); rm(m, idx); gc(verbose = FALSE)
}
counts <- new("dgCMatrix", i = i, p = p, x = x, Dim = c(length(genes), length(ord)), Dimnames = list(genes, ord))
rm(x, i, p); gc()
data <- counts; data@x <- dx; rm(dx); gc()

empty <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = dim(counts), dimnames = dimnames(counts))
so <- CreateSeuratObject(counts = empty, min.cells = 0, min.features = 0)
LayerData(so, layer = "data") <- empty
so@assays$RNA@layers$counts <- counts
so@assays$RNA@layers$data <- data
rm(counts, data, empty); gc()
validObject(so[["RNA"]])

md <- keep$meta[colnames(so), ]
md$orig.ident <- NULL
so@meta.data <- cbind(so@meta.data[, "orig.ident", drop = FALSE], md)
so$orig.ident <- factor(so$sample)
for (r in names(keep$red)) {                                # reorder embeddings to the object's cell order
  dr <- keep$red[[r]]; dr@cell.embeddings <- dr@cell.embeddings[colnames(so), , drop = FALSE]; so[[r]] <- dr
}
for (g in names(keep$graphs)) {
  gr <- keep$graphs[[g]]
  if (!identical(rownames(gr), colnames(so))) { gr <- as.Graph(gr[colnames(so), colnames(so)]); DefaultAssay(gr) <- "RNA" }
  so[[g]] <- gr
}
VariableFeatures(so) <- keep$hvg
Idents(so) <- "celltype"
so@misc <- keep$misc
so@misc$note <- paste(
  "GSE342497 (Edara et al. 2026, JCI Insight): infant rhesus macaque brain, ZIKV (n=2) vs Control (n=3), 14 dpi.",
  "Full transcriptome (Mmul10). counts = raw UMI; data = LogNormalize (scale 1e4). Harmony integration by sample.",
  "region = verified region (GEO sample titles swap Hippocampus/Striatum; file names are correct).",
  "Low-quality/doublet clusters removed. DE results for AQP4/DMD/DAG1 in misc$target_gene_DE.")
rm(keep); gc()
saveRDS(so, file.path(proj, "results/GSE342497_ZIKV_macaque_annotated_seurat.rds"))
message("cells: ", ncol(so), "  genes: ", nrow(so)); message("done")
