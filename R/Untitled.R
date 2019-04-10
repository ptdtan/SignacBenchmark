write_sp_sim <- function(f, d, prefix = "pbmc4k_batch2") {
  library(Matrix)
  sim <- readRDS(f)
  m <- as(sim$mat, "sparseMatrix")
  dir.create(file.path(paste0("data.2/", prefix)))
  dir.create(file.path(paste0("data.2/", prefix) , d))
  writeMM(m, file.path(paste0("data.2/", prefix), d, "matrix.mtx"))
  bx = make.unique(colnames(sim$mat))
  cluster <- data.frame(barcode = bx, cluster_id = c(rep(0, 500), rep(1, 500)))
  readr::write_tsv(cluster, file.path(paste0("data.2/", prefix), d, "cluster.tsv"))
  readr::write_tsv(data.frame(bx = bx), file.path(paste0("data.2/", prefix), d, "barcodes.tsv"), col_names = F)
  readr::write_tsv(data.frame(gene = rownames(sim$mat)), file.path(paste0("data.2/", prefix), d, "genes.tsv"), col_names = F)
}


load_de_result <- function(res, res_sim)
{
  hy <- readr::read_delim(res, delim = " ", col_names = F)[[1]][1:100]
  truth <- readRDS(res_sim)
  truth_genes <- c(truth[[1]]$low, truth[[2]]$low)
  inter.genes <- length(which(hy %in% truth_genes))
  res <- c(inter.genes/length(hy), inter.genes, hy)
  message("TPR: ", res[[1]], " inter.genes: ", inter.genes, " hy:", length(hy))
}

write_truth <- function(res_sim, dir)
{
  truth <- readRDS(res_sim)
  truth_genes <- truth_genes <- c(sapply(truth, function(s)s$low))
  readr::write_tsv(data.frame(gene = truth_genes), file.path(dir, "truth.tsv"), col_names = F)
}


load_stats_result <- function(f)
{
  stat <- readRDS(f)
  readr::write_tsv(as.data.frame(stat), paste0(f, ".tsv"))
}
