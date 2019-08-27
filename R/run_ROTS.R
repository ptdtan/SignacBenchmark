suppressPackageStartupMessages(library(ROTS))
suppressPackageStartupMessages(library(edgeR))

run_ROTScpm <- function(cells.1, cells.2) {
  message("ROTS, CPM")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      rots <- ROTS(data = mat, groups = as.numeric(as.factor(clusters)), B = 100, K = 500, log = FALSE, seed = 123)
    })

    list(session_info = session_info,
         timing = timing,
         res = rots,
         df = data.frame(pval = rots$pvalue,
                         padj = rots$FDR,
                         row.names = rownames(rots$data)))
  },  error = function(e) {
    "ROTS results could not be calculated"
    list(session_info = session_info,
         error = e)
  })
}
