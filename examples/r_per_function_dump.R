#!/usr/bin/env Rscript
# Per-function output dump for py-TSCAN's Notebook 3
# (function_by_function_R_parity.ipynb).
#
# Calls each public TSCAN R function in turn on lpsdata and saves its raw
# output to a per-function JSON file under examples/_r_outputs/. The Python
# notebook reads these and compares to its own calls.
#
# Run from py-TSCAN/ once before re-executing Notebook 3:
#   conda activate /scratch/users/steorra/env/CMAP
#   Rscript examples/r_per_function_dump.R

suppressMessages({
  library(TSCAN)
  library(jsonlite)
})

OUT_DIR <- "examples/_r_outputs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

data(lpsdata)
cat("[r-dump] lpsdata dims:", dim(lpsdata), "\n")

# ---- 1. preprocess --------------------------------------------------------
procdata <- preprocess(lpsdata,
                       takelog = TRUE,
                       logbase = 2,
                       pseudocount = 1,
                       minexpr_value = 1,
                       minexpr_percent = 0.5,
                       cvcutoff = 1)
write_json(
  list(
    shape = dim(procdata),
    rownames = rownames(procdata),
    colnames = colnames(procdata),
    head_3x3 = as.matrix(procdata)[1:3, 1:3],
    full     = as.matrix(procdata)
  ),
  file.path(OUT_DIR, "preprocess.json"),
  auto_unbox = TRUE, digits = NA, matrix = "rowmajor", na = "null", pretty = FALSE
)
cat("[r-dump] preprocess: dim =", dim(procdata), "\n")

# ---- 2. exprmclust --------------------------------------------------------
mc <- exprmclust(procdata,
                 clusternum = 2:9,
                 modelNames = "VVV",
                 reduce = TRUE)

# Also recompute the auto-selected pcadim — exprmclust hides it.
sdev <- prcomp(t(procdata), scale = TRUE)$sdev[1:20]
xv <- 1:20
optpoint <- which.min(sapply(2:10, function(i) {
  x2 <- pmax(0, xv - i)
  sum(lm(sdev ~ xv + x2)$residuals^2)
}))
pcadim <- optpoint + 1

# Capture MST edges (igraph → adjacency)
adj <- as.matrix(igraph::as_adjacency_matrix(mc$MSTtree))
edges <- which(upper.tri(adj) & adj > 0, arr.ind = TRUE)
mst_edges <- if (nrow(edges) > 0) cbind(edges[,1], edges[,2]) else matrix(integer(0), ncol=2)

write_json(
  list(
    pcadim = as.integer(pcadim),
    pcareduceres = as.matrix(mc$pcareduceres),
    pcareduceres_shape = dim(mc$pcareduceres),
    clusterid = as.integer(mc$clusterid),
    cluster_sizes = as.list(table(mc$clusterid)),
    clucenter = as.matrix(mc$clucenter),
    clucenter_shape = dim(mc$clucenter),
    mst_edges = mst_edges,
    cell_names = colnames(procdata)
  ),
  file.path(OUT_DIR, "exprmclust.json"),
  auto_unbox = TRUE, digits = NA, matrix = "rowmajor", na = "null", pretty = FALSE
)
cat("[r-dump] exprmclust: G =", length(unique(mc$clusterid)),
    ", pcadim =", pcadim, "\n")

# ---- 3. TSCANorder --------------------------------------------------------
ord_df <- TSCANorder(mc,
                     MSTorder = NULL,
                     orderonly = FALSE,
                     flip = FALSE,
                     listbranch = FALSE)

# Per-cell pseudotime keyed by procdata column order
all_cells <- colnames(procdata)
pt_per_cell <- rep(NA_real_, length(all_cells))
names(pt_per_cell) <- all_cells
pt_per_cell[ord_df$sample_name] <- ord_df$Pseudotime
state_per_cell <- rep(NA_integer_, length(all_cells))
state_per_cell[match(ord_df$sample_name, all_cells)] <- ord_df$State

write_json(
  list(
    ordered_names = as.character(ord_df$sample_name),
    State = as.integer(ord_df$State),
    Pseudotime = as.integer(ord_df$Pseudotime),
    n_ordered = nrow(ord_df),
    pseudotime_per_cell = as.numeric(pt_per_cell),
    state_per_cell = as.integer(state_per_cell)
  ),
  file.path(OUT_DIR, "TSCANorder.json"),
  auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE
)
cat("[r-dump] TSCANorder: ordered", nrow(ord_df), "cells\n")

# ---- 4. difftest (first 50 genes only — full takes minutes) --------------
de_subset <- procdata[1:50, ]
de <- difftest(de_subset, ord_df$sample_name, df = 3)
write_json(
  list(
    gene_names = rownames(de),
    pval = as.numeric(de$pval),
    FDR  = as.numeric(de$FDR),
    n_tested = nrow(de),
    sorted_order_genes = rownames(de)  # sort order R picks: (FDR, pval)
  ),
  file.path(OUT_DIR, "difftest.json"),
  auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE
)
cat("[r-dump] difftest: tested", nrow(de), "genes\n")

# ---- 5. orderscore --------------------------------------------------------
# Two-population coding: Unstimulated → 0, LPS_6h → 1
sub <- data.frame(
  cell = colnames(procdata),
  code = ifelse(grepl("Unstimulated", colnames(procdata)), 0L, 1L),
  stringsAsFactors = FALSE
)

# Compare 3 orderings: TSCAN order, reversed, sorted-by-name (deterministic "random")
random_order <- sample(sub$cell, length(sub$cell))  # use set.seed below
set.seed(42)
random_order <- sample(sub$cell, length(sub$cell))

orders <- list(ord_df$sample_name, rev(ord_df$sample_name), random_order)
scores <- orderscore(sub, orders)
write_json(
  list(
    scores = as.numeric(scores),
    labels = c("TSCAN", "reversed", "random_seed42"),
    sub_codes = sub$code,
    sub_cells = sub$cell,
    random_order = random_order
  ),
  file.path(OUT_DIR, "orderscore.json"),
  auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE
)
cat("[r-dump] orderscore: [", paste(round(scores, 4), collapse=", "), "]\n")

cat("\n[r-dump] all done. Files in", OUT_DIR, ":\n")
print(list.files(OUT_DIR))
