suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))

# run_edgeRLRT ------------------------------------------------------------

run_edgeRLRT <- function(cells.1, cells.2) {
  message("edgeRLRT")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      dge <- DGEList(mat, group = clusters)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~clusters)
      dge <- estimateDisp(dge, design = design)
      fit <- glmFit(dge, design = design)
      lrt <- glmLRT(fit)
      tt <- topTags(lrt, n = Inf)
    })
    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  },error = function(e) {
    "run_edgeRLRTrobust results could not be calculated"
    list(session_info = session_info,
         error = e)
  })
}

# run_edgeRQLFDetRate -----------------------------------------------------

run_edgeRQLFDetRate <- function(cells.1, cells.2) {
  message("edgeRQLFDetRate")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      dge <- DGEList(mat, group = clusters)
      dge <- edgeR::calcNormFactors(dge)
      cdr <- scale(colMeans(mat > 0))
      design <- model.matrix(~ cdr + clusters)
      dge <- estimateDisp(dge, design = design)
      fit <- glmQLFit(dge, design = design)
      qlf <- glmQLFTest(fit)
      tt <- topTags(qlf, n = Inf)
    })
    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  },error = function(e) {
    "run_edgeRLRTrobust results could not be calculated"
    list(session_info = session_info,
         error = e)
  })
}

# run_edgeRQLF ------------------------------------------------------------

run_edgeRQLF <- function(cells.1, cells.2) {
  message("edgeRQLF")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      dge <- DGEList(mat, group = clusters)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~clusters)
      dge <- estimateDisp(dge, design = design)
      fit <- glmQLFit(dge, design = design)
      qlf <- glmQLFTest(fit)
      tt <- topTags(qlf, n = Inf)
    })
    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  },error = function(e) {
    "run_edgeRLRTrobust results could not be calculated"
    list(session_info = session_info,
         error = e)
  })
}

# run_edgeRLRTrobust ------------------------------------------------------

run_edgeRLRTrobust <- function(cells.1, cells.2) {
  message("edgeRLRTrobust")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
      timing <- system.time({
      dge <- DGEList(mat, group = clusters)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~clusters)
      dge <- estimateGLMRobustDisp(dge, design = design)
      fit <- glmFit(dge, design = design)
      lrt <- glmLRT(fit)
      tt <- topTags(lrt, n = Inf)
    })

    plotBCV(dge)
    hist(tt$table$PValue, 50)
    hist(tt$table$FDR, 50)
    limma::plotMDS(dge, col = as.numeric(as.factor(clusters)), pch = 19)
    plotSmear(lrt)

    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  }, error = function(e) {
    "run_edgeRLRTrobust results could not be calculated"
    list(session_info = session_info,
         error = e)
  })
}


