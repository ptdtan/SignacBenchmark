run_Seurat <- function(cells.1, cells.2, method = "bimod") {
  message(paste("Seurat", method))
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      umicount <- mat.raw

      umicount <- umicount[, c(cells.1, cells.2)]
      colnames(umicount) <- c(rep("cond", length(cells.1)), rep("noncond", length(cells.2)))
      colnames(umicount) <- paste0(colnames(umicount), "__", 1:ncol(umicount))
      seur <- Seurat::CreateSeuratObject(raw.data = umicount, project = "scrnaseq", do.logNormalize = TRUE,
                                 names.field = 1, names.delim = "__")

      res <- Seurat::FindMarkers(seur, ident.1 = "cond",
                         ident.2 = "noncond", test.use = method, logfc.threshold = -1)

    })

    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$p_val,
                         padj = p.adjust(res$p_val, method = "BH"),
                         row.names = rownames(res)))
  }, error = function(e) {
    "Seurat results could not be calculated"
    list(session_info = session_info)
  })
}
