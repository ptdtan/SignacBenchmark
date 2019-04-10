suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))

run_MASTcpmDetRate <- function(cells.1, cells.2) {
  message("MAST, CPM (including detection rate)")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
    timing <- system.time({
      grp <- clusters
      cdr <- scale(colMeans(mat > 0))
      dge <- DGEList(counts = mat)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- edgeR::cpm(dge)
      sca <- FromMatrix(exprsArray = log2(cpms + 1),
                        cData = data.frame(wellKey = colnames(mat),
                                           grp = grp, cdr = cdr))
      zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
      mast <- lrTest(zlmdata, "grp")
    })

    list(session_info = session_info,
         timing = timing,
         res = mast,
         df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                         padj = mast[, "hurdle", "Pr(>Chisq)"],
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
  }, error = function(e) {
    "MASTcpmDetRate results could not be calculated"
    list(session_info = session_info)
  })
}

run_MASTcpm <- function(cells.1, cells.2) {
  message("MAST, CPM")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
    timing <- system.time({

      grp <- clusters
      dge <- DGEList(counts = mat)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- cpm(dge)
      sca <- FromMatrix(exprsArray = log2(cpms + 1),
                        cData = data.frame(wellKey = names(grp),
                                           grp = grp))
      zlmdata <- zlm.SingleCellAssay(~grp, sca)
      mast <- lrTest(zlmdata, "grp")
    })

    hist(mast[, "hurdle", "Pr(>Chisq)"], 50)

    list(session_info = session_info,
         timing = timing,
         res = mast,
         df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                         padj = mast[, "hurdle", "Pr(>Chisq)"],
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
  }, error = function(e) {
    "MASTcpmDetRate results could not be calculated"
    list(session_info = session_info)
  })

}

