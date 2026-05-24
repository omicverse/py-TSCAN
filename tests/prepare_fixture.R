#!/usr/bin/env Rscript
# One-time fixture preparation: dump TSCAN::lpsdata to RDS for both R and Python
# to consume. Run once from the CMAP conda env.

suppressMessages({
  library(TSCAN)
})

data(lpsdata)
cat("[fixture] lpsdata dims:", dim(lpsdata), "\n")
cat("[fixture] colname examples:", head(colnames(lpsdata), 3), "\n")

# Also dump as plain CSV so Python can load without R-RDS dependency.
out_rds <- "data/fixture_lpsdata.rds"
out_csv <- "data/fixture_lpsdata.csv"
saveRDS(lpsdata, out_rds)
write.csv(lpsdata, out_csv)
cat("[fixture] wrote:", out_rds, "and", out_csv, "\n")
