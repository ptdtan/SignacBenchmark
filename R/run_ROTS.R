suppressPackageStartupMessages(library(ROTS))
suppressPackageStartupMessages(library(edgeR))

run_ROTScpm <- function(cells.1, cells.2) {
  message("ROTS, CPM")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  timing <- system.time({
    grp <- clusters
    dge <- DGEList(counts = mat)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
    rots <- ROTS(data = cpms, groups = grp, B = 1000, K = 1000, log = FALSE, seed = 123)
  })

  hist(rots$pvalue, 50)
  hist(rots$FDR, 50)
  hist(rots$logfc, 50)
  print(rots$R)
  print(rots$Z)
  print(rots$k)
  print(rots$a1)
  print(rots$a2)

  list(session_info = session_info,
       timing = timing,
       res = rots,
       df = data.frame(pval = rots$pvalue,
                       padj = rots$pvalue,
                       row.names = rownames(rots$data)))
}
