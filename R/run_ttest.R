suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(genefilter))

run_ttest <- function(cells.1, cells.2) {
  message("t-test")

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  session_info <- sessionInfo()
  timing <- system.time({
    tmm <- edgeR::calcNormFactors(mat)
    tpmtmm <- edgeR::cpm(mat, lib.size = tmm * colSums(mat))
    logtpm <- log2(tpmtmm + 1)
    idx <- seq_len(nrow(logtpm))
    names(idx) <- rownames(logtpm)
    ttest_p <- sapply(idx, function(i) {
      t.test(logtpm[i, ] ~ clusters)$p.value
    })
  })

  list(session_info = session_info,
       timing = timing,
       df = data.frame(pval = ttest_p,
                       padj = ttest_p,
                       row.names = names(ttest_p)))
}
