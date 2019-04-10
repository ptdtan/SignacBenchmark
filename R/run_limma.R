suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

run_limmatrend <- function(cells.1, cells.2) {
  message("limmatrend")
  session_info <- sessionInfo()
  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      dge <- DGEList(mat, group = clusters)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~clusters)
      y <- new("EList")
      y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- lmFit(y, design = design)
      fit <- eBayes(fit, trend = TRUE, robust = TRUE)
      tt <- topTable(fit, n = Inf, adjust.method = "BH")
    })

    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$P.Value,
                         padj = tt$adj.P.Val,
                         row.names = rownames(tt)))
  },error = function(e) {
    "limmatrend results could not be calculated"
    list(session_info = session_info)
  })
}

run_voomlimma <- function(cells.1, cells.2) {
  message("voomlimma")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      dge <- DGEList(mat, group = clusters)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~clusters)
      vm <- voom(dge, design = design, plot = TRUE)
      fit <- lmFit(vm, design = design)
      fit <- eBayes(fit)
      tt <- topTable(fit, n = Inf, adjust.method = "BH")
    })

    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$P.Value,
                         padj = tt$adj.P.Val,
                         row.names = rownames(tt)))
  },error = function(e) {
    "voomlimma results could not be calculated"
    list(session_info = session_info)
  })
}
