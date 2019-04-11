process_islam <- function(root.dir)
{
  truth <- readr::read_tsv("Moliner_1000DEG.tsv")
  condition <- readr::read_csv("Moliner_condition.tsv")
  truth.genes <- truth$ensembl[!grepl("^NA", truth$ensembl)]
  truth.genes <- truth.genes[!grep("\.[1-9]$", truth.genes)]
  truth.genes <- truth.genes[!grepl("\\.[1-9]$", truth.genes)]
  dge <- readr::read_tsv("Islam/GSE29087.tsv")
  dge[[1]][is.na(dge[[1]])] <- "NA_"
  r.names <- dge[[1]]
  dge <- dge[, -1]
  dge <- as.matrix(dge)
  rownames(dge) <- r.names

  stopifnot(nrow(condition) == ncol(dge))
  data <- list(real = list(which(condition[[1]] == "Embryonic stem cell"),
                           which(condition[[1]] == "Embryonic fibroblast")),
               null = list(),
               mat = dge)
  truth <- list(ES_MEF = list(low = truth.genes,
                              high = c()))
  saveRDS(data, file.path(root.dir, paste0("ES_MEF.data.rds")))
  saveRDS(truth, file.path(root.dir, paste0("ES_MEF.truth.rds")))

  data <- list(real = list(which(condition[[1]] == "Embryonic fibroblast"),
                           which(condition[[1]] == "Embryonic stem cell")),
               null = list(),
               mat = dge)
  saveRDS(data, file.path(root.dir, paste0("MEF_ES.data.rds")))
  saveRDS(truth, file.path(root.dir, paste0("MEF_ES.truth.rds")))
}
