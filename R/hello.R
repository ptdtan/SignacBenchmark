
hello <- function() {
  print("Hello, world!")
}

library(Matrix)
library(Seurat)
#library(org.Hs.eg.db)


process_GSE115469 <- function()
{
  file <- "/Users/bioturing/Data/BenchMark/GSE115469/GSE115469_Data.csv"
  mat <- readr::read_csv(file)
  dir <- "/Users/bioturing/Data/BenchMark/GSE115469/"

  row.names <- mat[[1]]
  mat <- mat[, -1]
  col.names <- data.frame(colnames(mat))
  row.names.df <- select(org.Hs.eg.db, keys = row.names, keytype = "ALIAS", columns = "ENSEMBL")
  row.names <- row.names.df[match(row.names, row.names.df$ALIAS), ]
  genes_map <- readr::read_tsv("/Users/bioturing/Data/BenchMark/GSE115469/genes_map.tsv", col_names = F)
  row.names <- genes_map[[1]][match(row.names, genes_map$X2)]

  m <- as(mat, "matrix")
  m <- as(m, "dgCMatrix")
  rownames(m) <- row.names
  genes.df <- data.frame(id = row.names, symbol = genes_map[[2]][match(row.names, genes_map$X1)])

  for(i in seq(1, 5)){
    s_dir <- file.path(dir, paste0("P", i))
    dir.create(s_dir)
    cells.idx <- grep(paste0("^P", i), col.names[[1]])
    mm <- m[, cells.idx]
    barcodes <- data.frame(cell = col.names[[1]][cells.idx])
    readr::write_csv(barcodes, file.path(s_dir, "barcodes.tsv"), col_names = F)
    readr::write_tsv(genes.df, file.path(s_dir, "genes.tsv"), col_names = F)
    Matrix::writeMM(mm, file.path(s_dir, "matrix.mtx"))
  }

}

process_GSE62270 <- function()
{
  dir <- "/Users/bioturing/Data/BenchMark/GSE62270/"
  file.rds <- file.path(dir, "GSE62270-GPL17021.rds")

  download.file("http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE62270-GPL17021.rds",
                file.rds)
  sce <- readRDS(file.rds)
  sce.counts <- sce[[1]]
}

scoringGeneMarkers_cnt <- function(final_scores, zero.pct.max = 0.5,
                               max.genes = 30, plot = FALSE,
                               max.zero.cnt = 50,
                               min.case.pct = 0.5)
{
  final_scores <- final_scores[final_scores$`zero.pct` < zero.pct.max, ]

  final_scores$rate <- (final_scores$case.cnt)**2 * (final_scores$case.exp) *
    (1 / (final_scores$zero.pct + 0.001)) *
    (1 / (final_scores$zero.exp + 1))

  final_scores <- final_scores[order(-final_scores$rate), ]
  return(final_scores)
}

scoringGeneMarkers_rate <- function(final_scores, zero.pct.max = 0.5,
                                   max.genes = 30, plot = FALSE,
                                   max.zero.cnt = 50,
                                   min.case.pct = 0.5)
{
  #final_scores <- final_scores[final_scores$`zero.pct` < zero.pct.max, ]

  final_scores$rate <- (final_scores$case.cnt^2) * (final_scores$case.exp) *
    (1 / (final_scores$zero.pct + 0.001)) *
    (1 / (final_scores$zero.exp + 0.001))

  final_scores <- final_scores[order(-final_scores$rate), ]
  return(final_scores)
}
