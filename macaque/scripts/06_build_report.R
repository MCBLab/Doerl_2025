# Step 6: build the self-contained HTML report from the panel data written by 05.
# Everything is inlined: the statistics as a JSON literal, the supporting figures as base64 PNGs.
# The template (scripts/report_template.html) holds the prose and the chart code; this script only
# slims the data to what the charts read and substitutes the placeholders.
suppressPackageStartupMessages({ library(jsonlite); library(dplyr); library(magick) })
proj <- "/home/diego-coelho/Documents/mcblab/Edara2026"
d <- fromJSON(file.path(proj, "data/report_data.json"), simplifyDataFrame = TRUE)

# ---- columns of the heatmap: broad cell types plus the two astrocyte subtypes ----
cols <- c("Astrocytes (all)", "Astro GFAP-hi", "Astro proto", "Endothelial", "Pericytes",
          "Smooth muscle cells", "Microglia")
from_broad <- c("Astrocytes (all)" = "Astrocytes", "Endothelial" = "Endothelial", "Pericytes" = "Pericytes",
                "Smooth muscle cells" = "Smooth muscle cells", "Microglia" = "Microglia")
from_fine <- c("Astro GFAP-hi" = "Astrocytes_GFAP-hi", "Astro proto" = "Astrocytes_protoplasmic")

slim <- function(x, src, key, label) {
  x |> filter(.data[[key]] == src) |>
    transmute(gene, region, col = label, lfc = logFC, p = PValue, fdr = FDR,
              lcpm = log2CPM_ctrl, amb = ambient_prone, tested)
}
heat <- bind_rows(
  lapply(names(from_broad), function(l) slim(d$broad, from_broad[[l]], "celltype", l)),
  lapply(names(from_fine), function(l) slim(d$fine, from_fine[[l]], "celltype_fine", l)))
heat$lfc[!heat$tested] <- NA                       # a fold change without a test is not shown

groups <- lapply(names(d$panel), function(g) list(name = g, genes = d$panel[[g]]))

out <- list(
  cols = cols,
  groups = groups,
  heat = heat,
  cpm = d$cpm_astro,
  comp = d$comp |> select(animal, condition, region, celltype, pct),
  avpr = d$expr |> filter(gene == "AVPR1A") |> select(celltype, log2CPM, n_pb),
  built = format(Sys.Date(), "%e %B %Y") |> trimws())

# "</" inside a <script> would close it early; the escape is invisible to JSON.parse
json <- gsub("</", "<\\/", toJSON(out, auto_unbox = TRUE, na = "null", digits = 5), fixed = TRUE)

# ---- figures as base64 data URIs, downscaled to keep the file reasonable ----
img <- function(f, w = 1700) {
  p <- file.path(tempdir(), paste0("r_", basename(f)))
  image_write(image_resize(image_read(file.path(proj, "results/figures", f)), paste0(w, "x")), p, format = "png")
  paste0("data:image/png;base64,", base64_enc(readBin(p, "raw", file.size(p))))
}

html <- readLines(file.path(proj, "scripts/report_template.html"), warn = FALSE) |> paste(collapse = "\n")
html <- sub("/*__DATA__*/ null", json, html, fixed = TRUE)
html <- sub("__IMG_01__", img("01_UMAP_celltype_region_condition.png"), html, fixed = TRUE)
html <- sub("__IMG_05__", img("05_heatmap_log2FC_ZIKV_vs_Control.png", 1500), html, fixed = TRUE)
html <- sub("__IMG_07__", img("07_ambient_RNA_diagnostic.png", 1400), html, fixed = TRUE)
stopifnot(!grepl("__IMG_|__DATA__", html))

dest <- file.path(proj, "results/ZIKV_macaque_glymphatic_report.html")
writeLines(html, dest)
message("wrote ", dest, "  (", round(file.size(dest) / 1e6, 2), " MB)")
message("heat rows: ", nrow(heat), "  tested: ", sum(heat$tested))
