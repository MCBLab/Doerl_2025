# GSE342497 (Edara et al., JCI Insight 2026) — infant rhesus macaque brain, ZIKV vs control
# Step 1a: per-sample cell calling from raw matrices, QC, doublet removal, integration, clustering
suppressPackageStartupMessages({
  library(Seurat); library(DropletUtils); library(scDblFinder); library(harmony)
  library(dplyr); library(ggplot2); library(BiocParallel)
})
set.seed(42)
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
raw_dir <- file.path(proj, "data/raw")
out_dir <- file.path(proj, "results/qc"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 8e9)

files <- list.files(raw_dir, pattern = "raw_feature_bc_matrix.h5$", full.names = TRUE)
meta <- data.frame(file = files) |>
  mutate(base = basename(file),
         gsm = sub("_.*", "", base),
         animal = sub("^GSM[0-9]+_(M[0-9])_.*", "\\1", base),
         condition = ifelse(grepl("_Zika_", base), "ZIKV", "Control"),
         # region as encoded in the FILE NAME; GEO sample titles swap Hipp/Stri — verified with markers later
         region_file = unname(c(Amyg = "Amygdala", Hipp = "Hippocampus", Stri = "Striatum")[
           sub("^GSM[0-9]+_M[0-9]_[A-Za-z]+_([A-Za-z]+)_.*", "\\1", base)]),
         sample = paste(animal, condition, region_file, sep = "_"))
print(meta[, -1])

ckpt_dir <- file.path(proj, "data/per_sample_ckpt"); dir.create(ckpt_dir, showWarnings = FALSE)
# one sample per R process (keeps peak memory low): Rscript 01a_per_sample_qc.R <row index>
i <- as.integer(commandArgs(trailingOnly = TRUE)[1])
objs <- list(); qc_tab <- list()

  m <- meta[i, ]
  message("== ", m$sample)
  ckpt <- file.path(ckpt_dir, paste0(m$sample, ".rds"))
  if (file.exists(ckpt)) { message("checkpoint exists, skipping"); quit(save = "no") }
  raw <- Read10X_h5(m$file, use.names = TRUE)
  if (is.list(raw)) raw <- raw[["Gene Expression"]]
  # ZIKV samples were aligned to Mmul10 + ZIKV PRVABC59 (KU501215.1): 3 extra viral features.
  # Keep their per-barcode UMI total as metadata, then drop them so all samples share one gene list.
  viral <- intersect(c("upstream", "primer", "downstream"), rownames(raw))
  zikv_umi <- if (length(viral)) Matrix::colSums(raw[viral, , drop = FALSE]) else setNames(numeric(ncol(raw)), colnames(raw))
  raw <- raw[setdiff(rownames(raw), viral), ]
  # cell calling (raw/unfiltered matrices were deposited)
  ed <- emptyDrops(raw, lower = 200, niters = 10000, BPPARAM = SerialParam())
  keep <- which(ed$FDR <= 0.001)
  cnt <- raw[, keep]
  colnames(cnt) <- paste0(m$sample, "_", colnames(cnt))
  so <- CreateSeuratObject(cnt, project = m$sample, min.features = 0)
  so$sample <- m$sample; so$gsm <- m$gsm; so$animal <- m$animal
  so$condition <- m$condition; so$region_file <- m$region_file
  so$zikv_umi <- unname(zikv_umi[sub(paste0("^", m$sample, "_"), "", colnames(so))])
  so$percent.mt <- PercentageFeatureSet(so, features = intersect(
    rownames(so), c(grep("^MT-", rownames(so), value = TRUE),
                    # Mmul_10 mito genes are often unprefixed
                    c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB"))))
  so$percent.ribo <- PercentageFeatureSet(so, pattern = "^RP[SL][0-9]")
  n_called <- ncol(so)
  # adaptive QC: 3 MAD on log counts/features (lower), hard caps
  lc <- log10(so$nCount_RNA); lf <- log10(so$nFeature_RNA)
  pass <- so$nFeature_RNA >= 300 & lf > median(lf) - 3 * mad(lf) &
          lc < median(lc) + 4 * mad(lc) & so$percent.mt < 10
  so <- subset(so, cells = colnames(so)[pass])
  # doublets
  sce <- scDblFinder(as.SingleCellExperiment(so), BPPARAM = SerialParam())
  so$scDblFinder.class <- sce$scDblFinder.class; so$scDblFinder.score <- sce$scDblFinder.score
  qc_tab[[m$sample]] <- data.frame(sample = m$sample, called = n_called, pass_qc = ncol(so),
                                   doublets = sum(so$scDblFinder.class == "doublet"),
                                   median_umi = median(so$nCount_RNA), median_genes = median(so$nFeature_RNA),
                                   median_mt = median(so$percent.mt))
  objs[[m$sample]] <- subset(so, scDblFinder.class == "singlet")
  saveRDS(list(so = objs[[m$sample]], qc = qc_tab[[m$sample]]), ckpt)
message('sample done')
