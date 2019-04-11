run_all <- function(obj, type = "real", sample = F, n_samples = 100, seed = 1, res = list())
{
  source("R/run_Huy.R")
  source("R/run_Seurat.R")
  source("R/run_edgeR.R")
  source("R/run_MAST.R")
  source("R/run_limma.R")
  source("R/run_ttest.R")
  source("R/run_metagenomeseq.R")
  source("R/run_ROTS.R")
  source("R/run_SCDE.R")
  source("R/run_Wilcoxon.R")

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

  methods <- list(
                  res_SeuratT = function(...) run_Seurat(..., method = "t"),
                  res_SeuratPoisson = function(...) run_Seurat(..., method = "poisson"),
                  res_Seuratnegbinom = function(...) run_Seurat(..., method = "negbinom"),
                  res_edgeRLRT = run_edgeRLRT,
                  res_edgeRQLFDetRate = run_edgeRQLFDetRate,
                  res_edgeRQLF = run_edgeRQLF,
                  res_ttest = run_ttest,
                  res_MASTcpmDetRate = run_MASTcpmDetRate,
                  res_limmatrend = run_limmatrend,
                  run_Wilcoxon = run_Wilcoxon,
                  run_metagenomeSeq = run_metagenomeSeq,
                  run_MASTcpm = run_MASTcpm,
                  run_ROTScpm = run_ROTScpm,
                  run_voomlimma = run_voomlimma,
                  run_Harmony = run_Harmony
  )

  for (k in names(methods)) {
    if (is.null(res[[k]]))
      res[[k]] <- methods[[k]](cells.1, cells.2)
  }
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

run_TPR_real <- function(data_f, data.truth, prefix, sample = F, n_samples = 100, seed = 1)
{
  res.file <- file.path("real", paste(prefix, "res.rds", sep = "_"))
  truth <- readRDS(data.truth)
  obj = readRDS(data_f)
  s <- NULL
  res <- NULL
  tryCatch({
    timming <- system.time({
      if (file.exists(res.file)) {
        res.old = readRDS(res.file)
      } else{
        res.old = list()
      }
      res = run_all(obj, type = "real", sample = sample, n_samples = n_samples, seed = seed,
                    res = res.old)
      s = get_Precision_onedata(res, truth)
    })
    print(timming)
  }, error = function(e) {
    paste("Error", e)
  })

  dir.create("real")
  saveRDS(s, file = file.path("real", paste(prefix, "stats.rds", sep = "_")))
  saveRDS(res, file = file.path("real", paste(prefix, "res.rds", sep = "_")))
  return(s)
}

get_Precision_onedata <- function(res, truth)
{
  truth_genes <- c(sapply(truth, function(s)s$low))
  lapply(res, function(obj){
    genes <- rownames(obj$df)
    inter.genes <- length(which(genes[obj$df$padj <= 0.05] %in% truth_genes))
    res <- c(inter.genes/length(which(obj$df$padj <= 0.05)),
             inter.genes,
             length(which(obj$df$padj <= 0.05)))
  })
}

run_FDR <- function(file, prefix, type = "null")
{
    res.file <- file.path(type, paste0(prefix, ".stat.RDS"))
    obj = readRDS(file)
    fail = T
    tryCatch({
      timming <- system.time({
        if (file.exists(res.file)) {
          res.old = readRDS(res.file)
        } else{
          res.old = list()
        }
        res = run_all(obj, type = "null", res = res.old)
        s = get_FDR_onedata(res)
        saveRDS(res, res.file)
        fail = F
      })
      print(timming)
    }, error = function(e) {
      message("Run failed for", file, e)
      return(NULL)
    })
    if (fail)
      return(NULL)
    res <- cbind(s)
    colnames(res) <- prefix
    res <- data.frame(res[, 1], row.names = prefix)
    write.table(res, file.path(type, paste0(prefix, ".stats")))
}

draw_stats <- function(type = "null", output_folder = "figures",
                       title = "UMI", file.name = "FPR_UMI.svg")
{
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  files <- list.files(type)
  dfs <- lapply(file.path(type, files), read.table)
  rows <- lapply(dfs, unlist)
  stats <- do.call(rbind, rows)
  stats_m <- melt(stats)
  #stats_m[,3] <- stats_m[,3] + exp(-5)
  colnames(stats_m) <- c("x", "method", "FPR")
  new_methods <- factor(unlist(sapply(as.character(stats_m$method), function(m) {strsplit(m, "_")[[1]][2]})))
  stats_m$method <- factor(new_methods, levels = unique(new_methods))
  dftmp <- stats_m  %>%
    dplyr::filter(!is.na(FPR)) %>%
    dplyr::mutate(method = forcats::fct_reorder(method, FPR,
                                                fun = median, na.rm = FALSE,
                                                .desc = TRUE))
  cols = unique(dftmp$method)

  gglayers <- list(
    geom_hline(yintercept = 0.05),
    geom_boxplot(outlier.size = -1),
    geom_point(position = position_jitter(width = 0.2), size = 0.5),
    theme_bw(),
    xlab(""),
    ylab("FPR (fraction of genes with p < 0.05)"),
    scale_color_manual(values = cols),
    scale_y_sqrt(breaks = c(0.05, 0.5, 1), limits = c(0, 0.07)),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          legend.position = "none")
  )
  g <- ggplot(dftmp, aes(x = method, y = FPR, color = method)) + gglayers + ggtitle(title)

  ggsave(g, filename = paste0("figures/", file.name), width = 10, height = 5, units = "in")
  readr::write_tsv(stats_m, file.path(output_folder, paste0(type, "_stats.tsv") ))
}

