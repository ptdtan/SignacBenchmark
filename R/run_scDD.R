suppressPackageStartupMessages(library(scDD))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(SingleCellExperiment))


run_scDD <- function(cells.1, cells.2) {
  message("scDD")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
    timing <- system.time({
      sce <- SingleCellExperiment(assays=list(normcounts=mat),
                                  colData=data.frame(condition = clusters))
      prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
      scd <- scDD(sce, prior_param = prior_param, testZeroes = FALSE,
                  param = BiocParallel::MulticoreParam(workers = 2),
                  condition = "condition", min.size = 3, min.nonzero = NULL)
      res <- results(scd)
    })

    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$nonzero.pvalue,
                         padj = res$nonzero.pvalue.adj,
                         row.names = rownames(res)))
  }, error = function(e) {
    "scDD results could not be calculated"
    list(session_info = session_info)
  })
}
