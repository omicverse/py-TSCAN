#!/usr/bin/env Rscript
# Reference runner — invoked under conda env CMAP.
# Usage:  Rscript r_reference_driver.R <fixture.rds> <output.json>
#
# Produces the JSON the parity test reads. Keys MUST match
# `data/manifest.yaml::outputs[*].location_reference`.

suppressMessages({
  library(TSCAN)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
fixture_path <- args[1]
output_path  <- args[2]

if (is.na(fixture_path) || is.na(output_path)) {
  stop("Usage: Rscript r_reference_driver.R <fixture.rds> <output.json>")
}

# ---- load ------------------------------------------------------------------
lpsdata <- readRDS(fixture_path)
cat("[ref] fixture dims:", dim(lpsdata), "\n")

# ---- preprocess (exact defaults from TSCAN vignette) -----------------------
t0 <- proc.time()[["elapsed"]]
procdata <- preprocess(lpsdata)
t_pre <- proc.time()[["elapsed"]] - t0
cat("[ref] preprocess: ", t_pre, "s -> ", dim(procdata), " (genes x cells)\n")

# Also recompute the auto-selected pcadim for parity check (exprmclust hides it).
sdev <- prcomp(t(procdata), scale = TRUE)$sdev[1:20]
xv <- 1:20
optpoint <- which.min(sapply(2:10, function(i) {
  x2 <- pmax(0, xv - i)
  sum(lm(sdev ~ xv + x2)$residuals^2)
}))
pcadim <- optpoint + 1
cat("[ref] auto-selected pcadim =", pcadim, "\n")

# ---- exprmclust ------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
mc <- exprmclust(procdata)
t_mc <- proc.time()[["elapsed"]] - t0
cat("[ref] exprmclust: ", t_mc, "s; G =", length(unique(mc$clusterid)),
    "; pcareduceres dims:", dim(mc$pcareduceres), "\n")

# ---- TSCANorder ------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
ord_df <- TSCANorder(mc, orderonly = FALSE)
t_ord <- proc.time()[["elapsed"]] - t0
cat("[ref] TSCANorder: ", t_ord, "s -> ordered", nrow(ord_df), "cells\n")

# Build a per-cell pseudotime vector keyed by colnames(procdata),
# so Python can compare on the same cell axis.
all_cells <- colnames(procdata)
pseudotime_per_cell <- rep(NA_real_, length(all_cells))
names(pseudotime_per_cell) <- all_cells
pseudotime_per_cell[ord_df$sample_name] <- ord_df$Pseudotime

clusterid_per_cell <- as.integer(mc$clusterid[all_cells])

# ---- serialise to JSON -----------------------------------------------------
out <- list(
  cell_names        = all_cells,
  Pseudotime        = as.numeric(pseudotime_per_cell),
  State             = as.integer(ord_df$State[match(all_cells, ord_df$sample_name)]),
  ordered_names     = as.character(ord_df$sample_name),
  ordered_pseudotime = as.integer(ord_df$Pseudotime),
  clusterid         = clusterid_per_cell,
  pcadim            = as.integer(pcadim),
  pcareduceres      = as.matrix(mc$pcareduceres),
  clucenter         = as.matrix(mc$clucenter),
  procdata_dim      = as.integer(dim(procdata)),
  procdata_rownames = rownames(procdata),
  timings = list(
    preprocess = t_pre,
    exprmclust = t_mc,
    TSCANorder = t_ord
  )
)

write_json(out, output_path, auto_unbox = TRUE, digits = NA, matrix = "rowmajor",
           na = "null", pretty = FALSE)
cat("[ref] wrote", output_path, "\n")
