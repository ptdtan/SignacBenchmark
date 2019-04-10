suppressPackageStartupMessages(library(metagenomeSeq))

run_metagenomeSeq <- function(cells.1, cells.2) {
  message("metagenomeSeq")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
    timing <- system.time({
      obj <- newMRexperiment(mat,
                             phenoData = new("AnnotatedDataFrame",
                                             data = data.frame(condition = clusters,
                                                               row.names = colnames(mat))))
      p <- cumNormStatFast(obj)
      obj <- cumNorm(obj, p = p)
      mod <- model.matrix(~ condition, data = pData(obj))
      res <- fitFeatureModel(obj = obj, mod = mod, coef = 2)
      tbl <- MRtable(obj = res, number = Inf, by = 2, adjustMethod = "BH")
    })

    hist(tbl$pvalues, 50)

    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = tbl$pvalues,
                         padj = tbl$adjPvalues,
                         row.names = rownames(tbl)))
  }, error = function(e) {
    "metagenomeSeq results could not be calculated"
    list(session_info = session_info)
  })
}
