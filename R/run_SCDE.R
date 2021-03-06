suppressPackageStartupMessages(library(scde))

run_SCDE <- function(cells.1, cells.2) {
  message("scde")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  tryCatch({
    timing <- system.time({
    intcount <- apply(mat, 2, function(x) {storage.mode(x) <- 'integer'; x})
    o.ifm <- scde.error.models(counts = intcount, groups = clusters, n.cores = 56,
                               threshold.segmentation = TRUE,
                               save.crossfit.plots = FALSE, save.model.plots = FALSE,
                               verbose = 0, min.size.entries = min(300, nrow(mat) - 1))
    valid.cells <- o.ifm$corr.a > 0
    table(valid.cells)
    o.ifm <- o.ifm[valid.cells, ]
    o.prior <- scde.expression.prior(models = o.ifm, counts = intcount[, valid.cells],
                                     length.out = 400, show.plot = FALSE)
    grp <- factor(clusters[which(valid.cells)])
    names(grp) <- rownames(o.ifm)
    ediff <- scde.expression.difference(o.ifm, intcount[, valid.cells], o.prior,
                                        groups = grp, n.randomizations = 100,
                                        n.cores = 1, verbose = 0)
    p.values <- 2*pnorm(abs(ediff$Z), lower.tail = FALSE)
    p.values.adj <- 2*pnorm(abs(ediff$cZ), lower.tail = FALSE)
  })

  list(session_info = session_info,
       timing = timing,
       res = ediff,
       df = data.frame(pval = p.values,
                       padj = p.values.adj,
                       score = abs(ediff$Z),
                       row.names = rownames(ediff)))
  }, error = function(e) {
    "SCDE results could not be calculated"
    list(session_info = session_info)
  })
}
