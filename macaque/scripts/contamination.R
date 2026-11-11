# Shared ambient-RNA / doublet contamination estimate for pseudobulk changes (sourced by 03b and 05).
# ZIKV libraries carry more ambient RNA: e.g. ACTA2 rises 3-4 log2 in astrocytes, microglia and OPC (log2 CPM < 1)
# but only +0.4-0.7 in the smooth muscle cells that make it; astrocyte-only genes (AQP4) rise in non-astrocytes.
# For a tested group C (region x condition) the contamination ratio from source cell type S is the median over
# S-specific markers of CPM(marker in C) / CPM(marker in S); the expected foreign signal of gene g is
# E = sum_S ratio_S x CPM(g in S). A change is "ambient-prone" when it is not robust to removing E: the corrected
# log2FC, log2((O_zikv - E_zikv) / (O_ctrl - E_ctrl)) on pooled CPM O, keeps less than half of the uncorrected pooled
# log2FC (or flips sign), or E accounts for >= 90% of O in either condition. Markers can be genuinely expressed at
# low level in C, so E is an upper bound and the flag is conservative.
src_markers <- list(
  "Astrocytes" = c("GFAP", "SLC1A2", "SLC1A3", "GJA1", "ALDOC"),
  "Endothelial" = c("CLDN5", "FLT1", "PECAM1", "VWF", "ESAM"),
  "Smooth muscle cells" = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES"),
  "Pericytes" = c("RGS5", "KCNJ8", "ABCC9", "NOTCH3"),
  "Microglia" = c("P2RY12", "CX3CR1", "CSF1R", "C1QB", "TMEM119"),
  "Border-associated macrophages" = c("LYVE1", "MRC1", "F13A1"),
  "Fibroblasts" = c("DCN", "LUM", "COL1A2", "COL1A1"),
  "Oligodendrocytes" = c("MBP", "MOG", "PLP1", "MOBP"),
  "OPC" = c("PDGFRA", "CSPG4"),
  "Ependymal" = c("FOXJ1"))

# pooled CPM per group (columns named "region|condition|label"), restricted to `genes` + all source markers
pooled_cpm <- function(mat, grp, genes) {
  s <- sapply(split(seq_len(ncol(mat)), grp), function(j) rowSums(mat[, j, drop = FALSE]))
  (t(t(s) / colSums(s)) * 1e6)[union(genes, unlist(src_markers)), , drop = FALSE]
}

expected_foreign <- function(c_cpm, rg, cond, ct_broad, P_broad, genes) {
  expected <- 0
  for (S in setdiff(names(src_markers), ct_broad)) {
    key <- paste(rg, cond, S, sep = "|"); if (!key %in% colnames(P_broad)) next
    s_cpm <- P_broad[, key]; m <- src_markers[[S]][s_cpm[src_markers[[S]]] > 0]
    expected <- expected + median(c_cpm[m] / s_cpm[m]) * s_cpm[genes]
  }
  expected
}

# one row per gene for group (rg, ct); P holds ct's pooled CPM, P_broad the broad cell types (sources)
contamination <- function(genes, P, ct, ct_broad, rg, P_broad) {
  O <- sapply(c("Control", "ZIKV"), function(cd) P[genes, paste(rg, cd, ct, sep = "|")])
  E <- sapply(c("Control", "ZIKV"), function(cd)
    expected_foreign(P[, paste(rg, cd, ct, sep = "|")], rg, cd, ct_broad, P_broad, genes))
  O <- matrix(O, ncol = 2); E <- matrix(E, ncol = 2)
  lfc_obs <- log2((O[, 2] + 0.5) / (O[, 1] + 0.5))
  lfc_cor <- log2((pmax(O[, 2] - E[, 2], 0) + 0.5) / (pmax(O[, 1] - E[, 1], 0) + 0.5))
  data.frame(gene = genes, log2CPM_ctrl = log2(O[, 1] + 1),
             foreign_frac_ctrl = E[, 1] / O[, 1], foreign_frac_zikv = E[, 2] / O[, 2],
             log2FC_pooled = lfc_obs, log2FC_corrected = lfc_cor,
             ambient_prone = pmax(E[, 1] / O[, 1], E[, 2] / O[, 2]) >= 0.9 | lfc_cor / lfc_obs < 0.5)
}
