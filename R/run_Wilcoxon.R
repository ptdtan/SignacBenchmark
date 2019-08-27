suppressPackageStartupMessages(library(edgeR))

run_Wilcoxon <- function(cells.1, cells.2) {
  message("Wilcoxon")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
    timing <- system.time({
      tmm <- edgeR::calcNormFactors(mat)
      tpmtmm <- edgeR::cpm(mat, lib.size = tmm * colSums(mat))
      idx <- 1:nrow(tpmtmm)
      names(idx) <- rownames(tpmtmm)
      wilcox_p <- sapply(idx, function(i) {
        wilcox.test(tpmtmm[i, ] ~ clusters)$p.value
      })
    })

    hist(wilcox_p, 50)

    list(session_info = session_info,
         timing = timing,
         df = data.frame(pval = wilcox_p,
                         padj = p.adjust(wilcox_p),
                         row.names = names(wilcox_p)))
  }, error = function(e) {
    "Wilcoxon"
    list(session_info = session_info)
  })
}
