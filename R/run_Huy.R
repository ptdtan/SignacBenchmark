run_Huy <- function(cells.1, cells.2)
{
  message("Huy Marker")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time(
      {
        mat <- mat.raw[, c(cells.1, cells.2)]
        clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))
        scores <- constructGeneMarkersData(expr.mat = mat, celltypes = clusters, cell = "A",
                                          backend_Sum = DelayedMatrixStats::rowSums2,
                                          max.exp.zero = 0, case.thres = 0.2)
        res <- scoringGeneMarkers_rate(scores)
        rate <- log2(res$rate + 1)
        hist(rate,breaks = 50)
        m <- mean(rate)
        s <- sd(rate)
        res$pval <- 1 - pnorm(rate, mean = m, sd = s)
      })
    list(session_info = session_info,
         timing = timing,
         df = data.frame(pval = res$pval, padj = res$pval, row.names = res$name),
         res = res
         )
  }, error = function(e) {
    "Huy results could not be calculated"
    list(session_info = session_info)
  })
}

constructGeneMarkersData <- function (expr.mat, celltypes, cell, case.thres = 0.45, max.nonzero = 0.05,
          step = 10, zero.thres = 12, min.genes = 12, max.exp.zero = 1,
          backend_Sum = Matrix::rowSums)
{
  max.nonzero.cnt <- ncol(expr.mat) * max.nonzero
  is.cell <- celltypes == cell
  sum.cells <- sum(is.cell)
  is.cell <- which(is.cell == TRUE)
  nis.cell <- celltypes != cell
  sum.ncells <- sum(nis.cell)
  nis.cell <- which(nis.cell == TRUE)
  options(DelayedArray.block.size = 2e+09)
  t0 <- Sys.time()
  case.cnt <- backend_Sum(expr.mat[, is.cell, drop = F] >
                            0)
  good.genes <- which(case.cnt > case.thres * sum.cells)
  message(paste("Block row sums in ", Sys.time() - t0, "secs"))
  if (length(good.genes) == 0)
    return(NULL)
  case.cnt <- case.cnt[good.genes]
  A <- expr.mat[good.genes, , drop = F]
  t0 <- Sys.time()
  chunks <- round(nrow(A)/parallel::detectCores())
  zero.cnt <- backend_Sum(A[, nis.cell, drop = F] > max.exp.zero)
  case.exp <- log2(backend_Sum(A[, is.cell, drop = F])/sum(is.cell) + 1)
  zero.exp <- log2(backend_Sum(A[, nis.cell, drop = F])/sum(nis.cell) + 1)
  case.pct <- case.cnt/sum.cells
  zero.pct <- zero.cnt/sum.ncells
  total_genes <- rownames(expr.mat)[good.genes]
  scores.list <- list(case.pct = case.pct * 1, zero.pct = zero.pct *
                        1, zero.cnt = zero.cnt, case.exp = case.exp * 1, zero.exp = zero.exp,
                      case.cnt = case.cnt, name = total_genes)
  scores.df <- data.frame(scores.list)
}

constructGeneMarkersData.1 <- function (expr.mat, celltypes, cell, case.thres = 0.45, max.nonzero = 0.05,
                                      step = 10, zero.thres = 12, min.genes = 12, max.exp.zero = 1,
                                      backend_Sum = Matrix::rowSums)
{
  max.nonzero.cnt <- ncol(expr.mat) * max.nonzero
  is.cell <- celltypes == cell
  sum.cells <- sum(is.cell)
  is.cell <- which(is.cell == TRUE)
  nis.cell <- celltypes != cell
  sum.ncells <- sum(nis.cell)
  nis.cell <- which(nis.cell == TRUE)
  options(DelayedArray.block.size = 2e+09)
  t0 <- Sys.time()
  case.cnt <- backend_Sum(expr.mat[, is.cell, drop = F] >
                            0)
  good.genes <- which(case.cnt > case.thres * sum.cells)
  message(paste("Block row sums in ", Sys.time() - t0, "secs"))
  if (length(good.genes) == 0)
    return(NULL)
  case.cnt <- case.cnt[good.genes]
  A <- expr.mat[good.genes, , drop = F]
  t0 <- Sys.time()
  chunks <- round(nrow(A)/parallel::detectCores())
  zero.cnt <- backend_Sum(A[, nis.cell, drop = F] > max.exp.zero)
  case.exp <- backend_Sum(A[, is.cell, drop = F])/(case.cnt + 1)
  zero.exp <- backend_Sum(A[, nis.cell, drop = F])/(zero.cnt + 1)
  case.pct <- case.cnt/sum.cells
  zero.pct <- zero.cnt/sum.ncells
  total_genes <- rownames(expr.mat)[good.genes]
  scores.list <- list(case.pct = case.pct * 1, zero.pct = zero.pct *
                        1, zero.cnt = zero.cnt, case.exp = case.exp * 1, zero.exp = zero.exp,
                      case.cnt = case.cnt, name = total_genes)
  scores.df <- data.frame(scores.list)
}
