suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Matrix))

run_Harmony <- function(cells.1, cells.2) {
  message("Harmony")
  session_info <- sessionInfo()
  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  mat <- as(mat, "sparseMatrix")

  tryCatch({
    timing <- system.time({
    res <- Signac::HarmonyMarker(mat, as.numeric(as.factor(clusters)))
    })

    list(session_info = session_info,
         timing = timing,
         df = data.frame(pval = res$`P.adjusted.value`,
                         padj = res$`P.adjusted.value`,
                         row.names = res$`Gene.Name`))
  },error = function(e) {
    "Harmony results could not be calculated"
    list(session_info = session_info)
  })
}
