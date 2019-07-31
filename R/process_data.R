suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scDD))

process_null_pbmc4k <- function(root.dir, prefix, cluster = c(1,2,3,4,5))
{
  master_dir <- "data.1/PBMC4k"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = F)
  mat.raw <- as(mat.sp, "matrix")
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$graph$clusters

  for(i in cluster){
    message("Cluster", i)
    col.idx <- sample(which(clusters == i), 200, replace = F)
    mat_ <- mat.sp[, col.idx]

    null.data.1 <- seq(1, 100)
    null.data.2 <- seq(101, 200)
    PBMC4k <- list(real = list(),
                   null = list(null.data.1, null.data.2),
                   mat = as(mat_, "matrix"))
    data.file <- file.path(root.dir, paste0(prefix, "_null_", i, ".rds"))
    saveRDS(PBMC4k, data.file)
  }
}

process_null_macosko <- function(root.dir, prefix, cluster = c(1,2,3,4,5))
{
  master_dir <- "data.1/macosko2015/"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = F)
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$graph$clusters

  for(i in cluster){
    message("Cluster", i)
    col.idx <- sample(which(clusters == i), 200, replace = F)
    mat_ <- mat.sp[, col.idx]

    null.data.1 <- seq(1, 100)
    null.data.2 <- seq(101, 200)
    PBMC4k <- list(real = list(),
                   null = list(null.data.1, null.data.2),
                   mat = as(mat_, "matrix"))
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
                      nEE=4000, nEP=4000,
                      seed = 816)
{
  require(scDD)
  SD <- simulateSet(scDatEx, numSamples=numSamples, nDE=nDE, nDP=nDP, nDM=nDM,
                   nDB=nDB, nEE=nEE, nEP=nEP, sd.range=c(2,2), modeFC=c(2,3,4),
                   plots=FALSE,
                   random.seed=seed)
  return(SD)
}

  process_simulate <- function(output.folder, ndata = 3, nDE=500,
                               nDP=500, nDM=500, nDB=500,
                               nEE=9000, nEP=9000)
  {
    library(scDD)
    data(scDatEx)
    suppressPackageStartupMessages(library(SingleCellExperiment))
    suppressPackageStartupMessages(library(MultiAssayExperiment))
    for(i in seq(1, ndata)){
      file.data <- file.path(output.folder, paste0("sim_", i, "_data.rds"))
      file.truth <- file.path(output.folder, paste0("sim_", i, "_truth.rds"))
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

process_one_Fullength <- function(obj, meta.obj, root.dir, prefix, group, n.instances = 2)
{
  groups <- meta.obj[[group]]
  for (g in unique(groups)) {
    col.idx <- which(groups == unique(groups)[1])
    mat <- assay(obj)[, col.idx]
    generate_instance(mat, paste(prefix, g, sep="-"), root.dir, n.instances = n.instances)
  }
}

prepare_one_Fullength <- function(obj.file, meta.file, prefix = "a")
{
  obj <- readRDS(obj.file)
  gse <- GEOquery::getGEO(filename = meta.file)

  pdata <- Biobase::phenoData(gse)@data
  pdata <- pdata[, grep(":ch1", colnames(pdata))]
  sapply(colnames(pdata), function(col)table(pdata[[col]]))
  saveRDS(obj[[1]], paste0(prefix, ".rds"))
  saveRDS(pdata, paste0(prefix, "-metadata.rds"))
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
      prefix <- gsub(" ", "_", prefix)
      data.file <- file.path(dir, paste0(prefix, "_null_", s, ".rds"))
      saveRDS(obj, data.file)
  }
}

process_Conquer_UMI <- function()
{
  process_one_UMI("data/Conquer_UMI/10XMonoCytoT.rds", "data/Conquer_UMI/", "10XMonoCytoT", "group", 6)
  process_one_UMI("data/Conquer_UMI/UsoskinGSE59739.rds", "data/Conquer_UMI/", "UsoskinGSE59739", "Level.1", 3)
}


##### HOW TO GET FDR FOR NULL DATASETS.
#' All of the stuff here must be done in the results folder of conquer_comparison
#' Assume that we already have the results inside 'results' folder.
#' One .rds file for one method for a particular dataset
fdr_one_data <- function(d, method, filt = '') {
  if (filt==""){
	file.p <- paste0(paste(d, method, sep="_"), ".rds")
  } else{
	file.p <- paste0(paste(d, method, filt, sep="_"), ".rds")
  }
  if (!file.exists(file.p))
	  stop(paste("File", file.p, "doesn't exists, something went wrong!"))
  data <- readRDS(file.p)
  pval <- as.numeric(sapply(data, function(d){
    if (!"df" %in% names(d))
      return(NA)
    #col <- if (!"padj" %in% colnames(d)) "pval" else "padj"
    col <- "pval"
    s <- length(which(d$df[[col]] <= 0.05))
    return(s/(sum(!is.na(d$df[[col]])) + 0.001))}))
  return(pval)
}

#' Function to process one data
#'
#' @param methods_str string of methods, split by the comma
#' @param datasets list of datasets
#' @param file_name file name of the result RDS
#'
#' @return save file to the RDS file
#' @export
#'
#' @examples
process_one_set_data <- function(methods_str, datasets, file_name = "umi_non_filtering.RDS", filt = "") {
  #methods_str <- "BPSC,D3E,DESeq2betapFALSE,DESeq2census,DESeq2nofilt,DESeq2,DEsingle,edgeRLRTcensus,edgeRLRTdeconv,edgeRLRTrobust,edgeRLRT,edgeRQLFDetRate,edgeRQLF,limmatrend,MASTcpmDetRate,MASTcpm,MASTtpmDetRate,MASTtpm,metagenomeSeq,monoclecensus,monoclecount,monocle,NODESnofilt,NODES,ROTScpm,ROTStpm,ROTSvoom,SAMseq,scDD,SCDE,SeuratBimodIsExpr2,SeuratBimodnofilt,SeuratBimod,SeuratTobit,ttest,voomlimma,Wilcoxon"
  methods <- strsplit(methods_str, ",")[[1]]
  #datasets <- c("UsoskinGSE59739mock", "GSE62270-GPL17021mock", "10XMonoCytoTmock")
  pval_data_list <- list()
  for (m in methods){
    message("Method: ", m)
    pval_data_list[[m]] <- c()
    for (d in datasets) {
      message("Dataset: ", d)
      fdr <- fdr_one_data(d, m, filt)
      cat(paste(fdr))
      cat("\n")
      #message(m, " ", d, " ", paste(fdr))
      pval_data_list[[m]] <- c(pval_data_list[[m]], fdr)
    }
  }
  saveRDS(pval_data_list, file = file_name)
}


# methods_str <- "Venice,BPSC,D3E,DESeq2betapFALSE,DESeq2census,DESeq2nofilt,DESeq2,DEsingle,edgeRLRTcensus,edgeRLRTdeconv,edgeRLRTrobust,edgeRLRT,edgeRQLFDetRate,edgeRQLF,limmatrend,MASTcpmDetRate,MASTcpm,MASTtpmDetRate,MASTtpm,metagenomeSeq,monoclecensus,monoclecount,monocle,NODESnofilt,NODES,ROTScpm,ROTStpm,ROTSvoom,SAMseq,scDD,SCDE,SeuratBimodIsExpr2,SeuratBimodnofilt,SeuratBimod,SeuratTobit,ttest,voomlimma,Wilcoxon"
# #methods_str <-  "Venice"
# process_one_set_data(methods_str,
#                      c("UsoskinGSE59739mock", "GSE62270-GPL17021mock", "10XMonoCytoTmock"),
#                      "umi_non_filtering.RDS",
#                      "")
#
# process_one_set_data(methods_str,
#                      c("UsoskinGSE59739mock", "GSE62270-GPL17021mock", "10XMonoCytoTmock"),
#                      "umi_filtering_TPM_1_25p.RDS",
#                      "TPM_1_25p")
#
# process_one_set_data(methods_str,
#                     c("GSE45719mock", "GSE48968-GPL13112mock", "GSE60749-GPL13112mock", "GSE74596mock", "EMTAB2805mock"),
#                     "full-length_non_filtering.RDS",
#                     "")
#
# process_one_set_data(methods_str,
#                      c("GSE45719mock", "GSE48968-GPL13112mock", "GSE60749-GPL13112mock", "GSE74596mock", "EMTAB2805mock"),
#                      "full-length_filtering_TPM_1_25p.RDS",
#                      "")

#' Plot function for result from one RDS file
#'
#' @param rds.path
#' @param output_folder
#' @param title
#' @param file.name
#'
#' @return
#' @export
#'
#' @examples
plot_null_FPR_from_rds <- function(rds.path, output_folder = "figures",
                          title = "UMI non-filtering", file.name = "FPR_UMI.svg",
                          dataset = "UMI", filt = "non-filtering", return.data = FALSE,
                          total_methods_str = "BPSC,D3E,DESeq2betapFALSE,DESeq2census,DESeq2nofilt,DESeq2,DEsingle,edgeRLRTcensus,edgeRLRTdeconv,edgeRLRTrobust,edgeRLRT,edgeRQLFDetRate,edgeRQLF,limmatrend,MASTcpmDetRate,MASTcpm,MASTtpmDetRate,MASTtpm,metagenomeSeq,monoclecensus,monoclecount,monocle,NODESnofilt,NODES,ROTScpm,ROTStpm,ROTSvoom,SAMseq,scDD,SCDE,SeuratBimodIsExpr2,SeuratBimodnofilt,SeuratBimod,SeuratTobit,ttest,voomlimma,Wilcoxon")
{
  require(dplyr)
  require(ggplot2)
  #data <- readRDS("null_umi_non_filt.RDS")
  data <- readRDS(rds.path)
  for (d in names(data)) {
    if (all(is.na(data[[d]]))) {
      data[[d]] <- NULL
    } else{
      data[[d]] <- data[[d]][!is.na(data[[d]])]
    }
  }
  stats_m <- plyr::ldply(data, rbind)
  stats_m <- reshape2::melt(stats_m)
  stats_m <- stats_m[, c(1, 3)]
  colnames(stats_m) <- c("method", "FPR")

  stats_m[["dataset"]] <- dataset
  stats_m[["filt"]] <- filt
  return(stats_m)
}

plot_batch_null_data <- function()
{
  colors <- "#a65353,#4c0a00,#cc3600,#ffa280,#d9b1a3,#4d3e39,#995426,#402200,#f29d3d,#735c00,#f2ce3d,#59502d,#ccff00,#aab386,#739926,#39592d,#00f200,#00e600,#0d3312,#b6f2ce,#00ffaa,#269973,#435952,#006b73,#39dae6,#0d2b33,#3399cc,#7c98a6,#0081f2,#00294d,#0044ff,#203980,#bfd0ff,#00138c,#737899,#000033,#9979f2,#502d59,#302633,#e200f2,#f2b6ee,#b32d98,#ff0088,#73003d,#33001b,#806071,#d9003a,#ffbfc8"
  g_u_non_filter <- plot_null_FPR_from_rds(rds.path = "umi_non_filtering.RDS",
                                           dataset = "UMI",
                                           filt = "non-filtering",
                                           return.data = T) #+theme(axis.text.y = element_text(size = 9))
  g_f_non_filter <- plot_null_FPR_from_rds(rds.path = "full-length_non_filtering.RDS",
                                           dataset = "full-length",
                                           filt = "non-filtering",
                                           return.data = T)  #+theme(axis.text.y = element_text(size = 9))
  g_u_filter <- plot_null_FPR_from_rds(rds.path = "umi_filtering_TPM_1_25p.RDS",
                                       dataset = "UMI",
                                       filt = "filtering",
                                       return.data = T) #+theme(axis.text.y = element_text(size = 9))
  g_f_filter <- plot_null_FPR_from_rds(rds.path = "full-length_filtering_TPM_1_25p.RDS",
                                       dataset  = "full-length",
                                       filt     = "filtering",
                                       return.data = T)  #+theme(axis.text.y = element_text(size = 9))
  truefpr <- do.call(rbind, list(g_u_non_filter, g_f_non_filter, g_u_filter, g_f_filter))

  gglayers <- list(
    geom_hline(yintercept = 0.05),
    geom_boxplot(outlier.size = -1),
    geom_point(position = position_jitter(width = 0.2), size = 0.5),
    theme_bw(),
    xlab(""),
    ylab("FPR (fraction of genes with p < 0.05)"),
    scale_color_manual(values = strsplit(colors, ",")[[1]] ),
    scale_y_sqrt(breaks = c(0.05, 0.5, 1), limits = c(0, 1.25)),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          legend.position = "none")
          )

  g_non_filter <- truefpr %>% filter(filt=="non-filtering") %>%
    dplyr::mutate(method = forcats::fct_reorder(method, FPR,
                                                fun = median, na.rm = TRUE,
                                                .desc = TRUE)) %>%
    ggplot(aes(x = method, y = FPR, color = method)) +
    gglayers +
    stat_summary(fun.data = function(x) {
      return(data.frame(y = 1,
                        label = paste0("n=", sum(!is.na(x)))))},
      geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5,
      hjust = -0.2, angle = 90) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~ dataset, ncol = 1) + theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            strip.background = element_blank(),
                                            panel.border = element_rect(colour = "black"))

  g_filter <- truefpr %>% filter(filt=="filtering") %>%
    dplyr::mutate(method = forcats::fct_reorder(method, FPR,
                                                fun = median, na.rm = TRUE,
                                                .desc = TRUE)) %>%
    ggplot(aes(x = method, y = FPR, color = method)) +
    gglayers +
    stat_summary(fun.data = function(x) {
      return(data.frame(y = 1,
                        label = paste0("n=", sum(!is.na(x)))))},
      geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5,
      hjust = -0.2, angle = 90) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~ dataset, ncol = 1) + theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            strip.background = element_blank(),
                                            panel.border = element_rect(colour = "black"))

  gg <- cowplot::plot_grid(g_non_filter + ggtitle("Without filtering"),
                           g_filter + ggtitle("After filtering"))
  pdf("Null_figure1.pdf", width = 12, height = 7)
  print(gg)
  dev.off()

  ggsave(gg, width = 12, height = 7, filename = "Null_figure1.png")
}
