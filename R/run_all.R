run_all <- function(obj, type = "real", sample = F, n_samples = 100, seed = 1)
{
  source("R/run_Huy.R")
  source("R/run_Seurat.R")
  source("R/run_edgeR.R")
  source("R/run_MAST.R")
  source("R/run_limma.R")
  source("R/run_ttest.R")
  mat.raw <<- obj[["mat"]]
  set.seed(seed)
  cells.1 <- obj[[type]][[1]]
  cells.2 <- obj[[type]][[2]]

  n.samples.1 <- sample(seq(25, 50))
  if(sample){
    cells.1 <- sample(cells.1, min(length(cells.1), n.samples.1), replace = F)
    cells.2 <- sample(cells.2, min(length(cells.2), 1000), replace = F)
  }

  message("cells.1 ", length(cells.1))
  message("cells.2 ", length(cells.2))

  stopifnot(length(intersect(cells.1, cells.2)) == 0)
  res <- list(res_Huy = run_Huy(cells.1, cells.2),
              res_SeuratBimod = run_Seurat(cells.1, cells.2, method = "bimod"),
              res_SeuratT = run_Seurat(cells.1, cells.2, method = "t"),
              res_SeuratPoisson = run_Seurat(cells.1, cells.2, method = "poisson"),
              res_Seuratnegbinom = run_Seurat(cells.1, cells.2, method = "negbinom"),
              res_edgeRLRT = run_edgeRLRT(cells.1, cells.2),
              res_edgeRQLFDetRate = run_edgeRQLFDetRate(cells.1, cells.2),
              res_edgeRQLF = run_edgeRQLF(cells.1, cells.2),
              res_ttest = run_ttest(cells.1, cells.2),
              res_MASTcpmDetRate = run_MASTcpmDetRate(cells.1, cells.2),
              res_limmatrend = run_limmatrend(cells.1, cells.2),
              run_Wilcoxon = run_Wilcoxon(cells.1, cells.2)
              )
  return(res)
}

run_FDR_sim <- function()
{
  list_files <- paste0("data/sim_", seq(1, 10),".rds")
  obj <- readRDS(list_files[1])
}

run_FDR_real <- function()
{
  files <- list(
    PBMC4k_1_sim = "data.1/PBMC4k_1_sim_data.rds"
    # GSE62270 = "data.1/GSE62270_data.rds",
    # GSE81076_GPL16791 = "data.1/GSE81076_GPL16791_data.rds",
    # GSE81076_GPL18573 = "data.1/GSE81076_GPL18573_data.rds",
    # pbmc4k_1 = "data.1/PBMC4k_1_data.rds",
    # pbmc4k_2 = "data.1/PBMC4k_2_data.rds"
    # pbmc4k_3 = "data.1/PBMC4k_3_data.rds",
    # zeisel2015_7 = "data.1/zeisel2015_7_data.rds",
    # zeisel2015_8 = "data.1/zeisel2015_8_data.rds",
    # zeisel2015_9 = "data.1/zeisel2015_9_data.rds"
    )
  stats <- lapply(files, function(file){
                  obj = readRDS(file)
                  tryCatch({
                    timming <- system.time({
                      res = run_all(obj, type = "real", sample = F)
                      s = get_FDR_onedata(res)
                    })
                    print(timming)
                    return(s)
                }, error = function(e) {
                  paste("Run failed for", file)
                  })
              })
  saveRDS(stats, file = "data.1/stats_results.rds")
  return(stats)
}

get_FDR_onedata <- function(res)
{
  lapply(res, function(obj){
    length(which(obj$df$padj < 0.05))/nrow(obj$df)
  })
}

run_TPR_real <- function(data_f, data.truth, prefix, sample = F, n_samples = 100, seed =1)
{
  truth <- readRDS(data.truth)
  obj = readRDS(data_f)
  s <- NULL
  res <- NULL
  tryCatch({
    timming <- system.time({
      res = run_all(obj, type = "real", sample = sample, n_samples = n_samples, seed = seed)
      s = get_Precision_onedata(res, truth)
    })
    print(timming)
  }, error = function(e) {
    paste("Error", e)
  })

  saveRDS(s, file = paste(data_f, prefix, "stats.rds", sep = "_"))
  saveRDS(res, file = paste(data_f, prefix, "res.rds", sep = "_"))
  return(s)
}

get_Precision_onedata <- function(res, truth)
{
  truth_genes <- c(sapply(truth, function(s)s$low))
  lapply(res, function(obj){
    genes <- rownames(obj$df)
    inter.genes <- length(which(genes[obj$df$padj <= 0.05] %in% truth_genes))
    res <- c(inter.genes/length(which(obj$df$padj <= 0.05)), inter.genes, length(which(obj$df$padj <= 0.05)))
  })
}

run_one_study <- function(file, type, prefix)
{
    obj = readRDS(file)
    tryCatch({
      timming <- system.time({
        res = run_all(obj, type = type)
        s = get_FDR_onedata(res)
      })
      print(timming)
    }, error = function(e) {
      message("Run failed for", file)
      return(NULL)
    })
    res <- cbind(s)
    colnames(res) <- prefix
    res <- data.frame(res[, 1], row.names = prefix)
    write.table(res, file.path(type, paste0(prefix, ".stats")))
}

draw_stats <- function(type = "null", output_folder = "figures")
{
  require(ggplot2)
  require(reshape2)
  files <- list.files(type)
  dfs <- lapply(file.path(type, files), read.table)
  rows <- lapply(dfs, unlist)
  stats <- do.call(rbind, rows)
  stats_m <- melt(stats)
  readr::write_tsv(stats_m, file.path(output_folder, paste0(type, "_stats.tsv") ))
}
