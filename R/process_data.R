suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scDD))

process_pbmc4k <- function(root.dir, prefix, cluster = c(1,2,3))
{
  master_dir <- "data.1/PBMC4k"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  mat.raw <- as(mat.sp, "matrix")
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$kmeans$clusters[3, ]

  for(i in cluster){
    message("Cluster", i)
    null.data <- which(clusters == i)
    null.data.1 <- sample(x = null.data, min(100, length(null.data)), replace = F)
    d <- setdiff(null.data, null.data.1)
    null.data.2 <- sample(x = d ,
                          min(100, length(d)), replace = F)
    real.data.1 <- which(clusters == i)
    real.data.2 <- which(clusters != i)
    PBMC4k <- list(real = list(real.data.1, real.data.2),
                     null = list(null.data.1, null.data.2),
                     mat = assay(mat.raw))
    data.file <- file.path(root.dir, paste0(prefix, "_null_", i, ".rds"))
    saveRDS(PBMC4k, data.file)
  }
}

process_zeisel2015 <- function()
{
  master_dir <- "data.1/zeisel2015"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  mat.raw <- as(mat.sp, "matrix")
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$graph$clusters
  cluster <- c(7,8,9)

  for(i in cluster){
    message("Cluster", i)
    null.data <- which(clusters == i)
    null.data.1 <- sample(x = null.data, min(100, sum(clusters == i)), replace = F)
    d <- setdiff(null.data, null.data.1)
    null.data.2 <- sample(x = d ,
                          min(100, length(d)), replace = F)
    real.data.1 <- which(clusters == i)
    real.data.2 <- which(clusters != i)
    data <- list(real = list(real.data.1, real.data.2),
                   null = list(null.data.1, null.data.2),
                   mat = assay(mat.raw))
    data.file <- file.path(paste0("data.1/zeisel2015_", i, "_data.rds"))
    saveRDS(data, data.file)
  }
}

process_GSE62270 <- function()
{
  file <- "data.1/GSE62270-GPL17021.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "source_name_ch1",
    keepgroups = c("Randomly extracted cells from whole intestinal organoids",
                   "Randomly extracted ex vivo isolated 5 day YFP positive cells")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "Randomly extracted Lgr5-positive intestinal cells"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  data.file <- file.path("data.1/GSE62270_data.rds")

  GSE62270 <- list(real = list(real.data.1, real.data.2),
                   null = list(null.data.1, null.data.2),
                   mat = assay(data[[1]]))
  saveRDS(GSE62270, data.file)
}

process_GSE81076 <- function()
{

  # GSE81076_GPL18573 -------------------------------------------------------
  file <- "data.1/GSE81076-GPL18573.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "characteristics_ch1",
    keepgroups = c("cell type: CD13+ sorted cells",
                   "cell type: CD24+ CD44+ live sorted cells")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "cell type: live sorted cells"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  GSE81076_GPL18573 <- list(real = list(real.data.1, real.data.2),
                  null = list(null.data.1, null.data.2),
                  mat = assay(data[[1]]))
  data.file <- file.path("data.1/GSE81076_GPL18573_data.rds")
  saveRDS(GSE81076_GPL18573, data.file)

  # GSE81076-GPL16791.rds ---------------------------------------------------
  file <- "data.1/GSE81076-GPL16791.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "characteristics_ch1",
    keepgroups = c("cell type: exocrine fraction, live sorted cells",
                   "cell type: live sorted cells")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "cell type: live sorted cells"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  GSE81076_GPL16791 <- list(real = list(real.data.1, real.data.2),
                            null = list(null.data.1, null.data.2),
                            mat = assay(data[[1]]))
  data.file <- file.path("data.1/GSE81076_GPL16791_data.rds")
  saveRDS(GSE81076_GPL16791, data.file)

}

process_GSE78779 <- function()
{
  file <- "data.1/GSE78779-GPL17021.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "characteristics_ch1.1",
    keepgroups = c("tissue: Ear",
                   "tissue: femora and tibiae")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "tissue: Ear"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  GSE78779 <- list(real = list(real.data.1, real.data.2),
                   null = list(null.data.1, null.data.2),
                   mat = assay(data[[1]]))
  data.file <- file.path("data.1/GSE78779_data.rds")
  saveRDS(GSE78779, data.file)
}

#' Title
#'
#' @param nDE
#' @param nDP
#' @param nDM
#' @param nDB
#' @param nEE
#' @param nEP
#' @param numSamples
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
simulate <-  function(scDatEx, numSamples=500,
                      nDE=500, nDP=500, nDM=500, nDB=500,
                      nEE=9000, nEP=9000,
                      seed = 816)
{
  SD <- simulateSet(scDatEx, numSamples=numSamples, nDE=nDE, nDP=nDP, nDM=nDM,
                   nDB=nDB, nEE=nEE, nEP=nEP, sd.range=c(2,2), modeFC=c(2,3,4),
                   plots=FALSE,
                   random.seed=seed)
  return(SD)
}

process_simulate <- function()
{
  library(scDD)
  data(scDatEx)
  for(i in seq(1, 3)){
    file.data <- file.path(paste0("data.2/sim_", i, "_data.rds"))
    file.truth <- file.path(paste0("data.2/sim_", i, "_truth.rds"))
    message(paste("Doing", file.data))
    SD <- simulate(scDatEx)
    real.data <- colData(SD)[["condition"]]
    real.data.1 <- which(real.data == 1)
    real.data.2 <- which(real.data != 1)
    sim <- list(real = list(real.data.1, real.data.2),
                null = NULL,
                mat = assay(SD))
    saveRDS(sim, file = file.data)
    truth <- list(total = list(low = rownames(sim$mat)[grep("^D", rownames(sim$mat))]))
    saveRDS(truth, file = file.truth)
  }
}

prepare_all <- function()
{
  message("PBMC4k")
  process_pbmc4k()

  message("process_zeisel2015")
  process_zeisel2015()

  message("process_GSE62270")
  process_GSE62270()

  message("process_GSE81076")
  process_GSE81076()

  message("process_GSE78779")
  process_GSE78779()
}

install_packages <- function()
{
  packages <- c("edgeR", "MultiAssayExperiment", "scDD", "genefilter")
  BiocInstaller::biocLite(packages)
}


process_one_UMI <- function(file.obj, root.dir, prefix, group, n.instances = 3)
{
  obj <- readRDS(file.obj)
  obj <- BiocGenerics:::updateObject(obj)
  groups <- colData(obj)[[group]]
  col.idx <- which(groups == unique(groups)[1])
  mat <- assay(obj)[, col.idx]
  generate_instance(mat, prefix, root.dir, n.instances = n.instances)
}

generate_instance <- function(mat, prefix, dir, n.instances = 3) {
  dir.create(dir)
  samples.sizes <- c()
  for (i in seq(0, n.instances - 1)) {
    samples.sizes <- c(samples.sizes, as.integer(ncol(mat)/(i + 2)))
  }

  for (s in samples.sizes) {
      message(dir, " ", prefix, " ", s)
      select.idx <- sample(seq(1, ncol(mat)), s, replace = F)
      obj <- list(real = NULL,
                  null = list(select.idx, setdiff(seq(1, ncol(mat)), (select.idx))),
                  mat = mat)
      data.file <- file.path(dir, paste0(prefix, "_null_", s, ".rds"))
      saveRDS(obj, data.file)
  }
}

process_Conquer_UMI <- function()
{
  process_one_UMI("data/Conquer_UMI/10XMonoCytoT.rds", "data/Conquer_UMI/", "10XMonoCytoT", "group", 6)
  process_one_UMI("data/Conquer_UMI/UsoskinGSE59739.rds", "data/Conquer_UMI/", "UsoskinGSE59739", "Level.1", 3)
}
