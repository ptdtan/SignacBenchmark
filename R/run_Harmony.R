suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Matrix))

run_Harmony <- function(cells.1, cells.2) {
  message("Venice")
  session_info <- sessionInfo()
  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  mat <- as(mat, "sparseMatrix")

  tryCatch({
    timing <- system.time({
      dge <- DGEList(counts = mat)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- edgeR::cpm(dge)
    })

  res <- Signac::VeniceMarker(mat, as.numeric(as.factor(clusters)) - 1)
  list(session_info = session_info,
       timing = timing,
       df = data.frame(pval = 10^res$`Log10.p.value`,
                       padj = 10^res$`Log10.adjusted.p.value`,
                       row.names = res$`Gene.Name`))
  },error = function(e) {
    "Venice results could not be calculated"
    list(session_info = session_info)
  })
}
