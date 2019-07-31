run_all <- function(obj, type = "real", sample = F, seed = 1, res = list(),
                    filtering = TRUE)
{
  # source("R/run_Huy.R")
  # source("R/run_Seurat.R")
  # source("R/run_edgeR.R")
  # source("R/run_MAST.R")
  # source("R/run_limma.R")
  # source("R/run_ttest.R")
  # source("R/run_metagenomeseq.R")
  # source("R/run_ROTS.R")
  # source("R/run_SCDE.R")
  # source("R/run_Wilcoxon.R")

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

  ### Filtering
  if (filtering) {
    message("Filtering")
    rowsum.cell12 <- DelayedMatrixStats::rowSums2(mat.raw[, c(cells.1, cells.2)] > 0)
    mat.raw <- mat.raw[which(rowsum.cell12 > 0.25 * (length(c(cells.1, cells.2)))), ]
  }
  methods <- list(
                  res_SeuratT = function(...) run_Seurat(..., method = "t"),
                  res_SeuratPoisson = function(...) run_Seurat(..., method = "poisson"),
                  res_Seuratnegbinom = function(...) run_Seurat(..., method = "negbinom"),
                  res_SeuratWilcox = function(...) run_Seurat(..., method = "wilcox"),
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

run_TPR_real <- function(data_f, data.truth, prefix, sample = F,
                          seed = 1, filtering = TRUE)
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
      res = run_all(obj, type = "real", sample = sample, seed = seed,
                    res = res.old, filtering = )
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
    idx <- which(obj$df$padj <= 0.05)
    inter.genes <- length(which(genes[idx] %in% truth_genes))
    res <- c(inter.genes/length(idx),
             inter.genes,
             length(idx))
  })
}

run_FDR <- function(file, prefix, out.dir, type = "null", filtering = TRUE)
{
    dir.create(out.dir)
    res.file <- file.path(out.dir, paste0(prefix, ".stat.RDS"))
    obj = readRDS(file)
    fail = T
    tryCatch({
      timming <- system.time({
        if (file.exists(res.file)) {
          res.old = readRDS(res.file)
        } else{
          res.old = list()
        }
        res = run_all(obj, type = "null", res = res.old, filtering = filtering)
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
    write.table(res, file.path(out.dir, paste0(prefix, ".stats")))
}

plot_null_FPR <- function(input_folder = "null", output_folder = "figures",
                       title = "UMI", file.name = "FPR_UMI.svg")
{
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  files <- list.files(input_folder)
  files <- files[grep("stats", files)]
  dfs <- lapply(file.path(input_folder, files), read.table)
  rows <- lapply(dfs, unlist)
  stats <- do.call(rbind, rows)
  stats_m <- melt(stats)
  #stats_m[,3] <- stats_m[,3] + exp(-5)
  colnames(stats_m) <- c("x", "method", "FPR")
  new_methods <- unlist(sapply(as.character(stats_m$method), function(m) {strsplit(m, "_")[[1]][2]}))
  stats_m$method <- factor(new_methods, levels = unique(new_methods))
  dftmp <- stats_m  %>%
    dplyr::filter(!is.na(FPR)) %>%
    dplyr::mutate(method = forcats::fct_reorder(method, FPR,
                                                fun = mean, na.rm = FALSE,
                                                .desc = TRUE))
  cols = unique(dftmp$method)

  colors_str <- '#a65353, #4c0a00, #cc3600, #ffa280, #d9b1a3, #4d3e39, #995426, #402200, #f29d3d, #735c00, #f2ce3d, #59502d, #ccff00, #aab386, #739926, #39592d, #00f200, #00e600, #0d3312, #b6f2ce, #00ffaa, #269973, #435952, #006b73, #39dae6, #0d2b33, #3399cc, #7c98a6, #0081f2, #00294d, #0044ff, #203980, #bfd0ff, #00138c, #737899, #000033, #9979f2, #502d59, #302633, #e200f2, #f2b6ee, #b32d98, #ff0088, #73003d, #33001b, #806071, #d9003a, #ffbfc8'
  colors <- strsplit(colors_str, ', ')[[1]]

  gglayers <- list(
    geom_hline(yintercept = 0.05),
    geom_boxplot(outlier.size = 1),
    geom_point(position = position_jitter(width = 0.2), size = 0.5),
    theme_bw(),
    xlab(""),
    ylab(""),
    scale_fill_manual(values = colors ),
    scale_y_sqrt(limits = c(0, 1)),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 11),
          legend.position = "none")
  )
  g <- ggplot(dftmp, aes(x = method, y = FPR, fill = method)) + gglayers + ggtitle(title)

  ggsave(g, filename = paste0("figures/", file.name), width = 10, height = 5, units = "in")
  readr::write_tsv(stats_m, file.path(output_folder, paste0("null", "_stats.tsv") ))
  return(g)
}

plot_null_4fig <- function()
{
  g_u_non_filter <- plot_null_FPR(input_folder = "null/UMI/", output_folder = "figures",
                                  title = "UMI", file.name = "FPR_UMI.png") + ylab("FPR (fraction of genes with p < 0.05)") +theme(axis.text.y = element_text(size = 9))
  g_f_non_filter <- plot_null_FPR(input_folder = "null/full-length/", output_folder = "figures",
                                  title = "full-length", file.name = "FPR_full-length.png")
  g_u_filter <- plot_null_FPR(input_folder = "null/UMI-filtering/", output_folder = "figures",
                                  title = "UMI-filtering", file.name = "FPR_UMI-filtering.png") + ylab("FPR (fraction of genes with p < 0.05)") +theme(axis.text.y = element_text(size = 9))

  g_f_filter <- plot_null_FPR(input_folder = "null/full-length-filtering/", output_folder = "figures",
                              title = "full-length-filtering", file.name = "FPR_full-length-filtering.png")
  gg <- cowplot::plot_grid(g_u_non_filter, g_f_non_filter, g_u_filter, g_f_filter)
  ggsave(gg, width = 10, height = 7, filename = "FPR_NULL.svg", units = "in")
}

plot_ROC_AUC <- function(res_f, truth_f, prefix = "ES_MEF_ROC_cpm.svg", xlim = c(0, 1))
{
  res <- readRDS(res_f)
  truth <- readRDS(truth_f)

  truth.genes <- c(sapply(truth, function(s)s$low))
  library(pROC)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  sensitivities <- c()
  one_m_specification <- c()
  methods <- c()
  auc <- c()
  for (i in seq(1, length(res))){
    test <- res[[i]]$df
    pred <- 1 - test$padj
    y <- match(rownames(test), truth.genes)
    y <- as.numeric(is.na(y))
    r <- roc(y, pred, direction = "auto")
    sensitivities <- c(sensitivities, r$sensitivities)
    one_m_specification <- c(one_m_specification, 1 - r$specificities)
    m <- names(res)[i]
    m <- strsplit(m, "_")[[1]][2]
    methods <- c(methods, rep(paste0(m, " AUC=", round(r$auc, 3)),
                              length(r$sensitivities)))
    auc <- c(auc, rep(r$auc, length(r$sensitivities)))
  }
  data <- data.frame(sens = sensitivities,
                     one_m_spec = one_m_specification,
                     method = methods,
                     auc = auc) %>% arrange(-auc)
  data$method <- factor(data$method, levels = unique(data$method))
  colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000')
  g <- ggplot(data,aes(one_m_spec, sens, color=method))+geom_line(size = 0.8, alpha = 0.7)+
    labs(title= "ROC curve",
         x = "False Positive Rate (1-Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    scale_x_continuous(limits = xlim) +
    theme(legend.position = c(0.67, 0.18),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10))+
    guides(col = guide_legend(ncol = 2))
  ggsave(g, filename = file.path("figures", prefix), width = 7, height = 5, units = "in", dpi = 120, device = "svg")
}

plot_frac_DE_genes <- function(res_f, truth_f)
{
  require(ggplot2)
  require(dplyr)
  res <- readRDS(res_f)
  truth <- readRDS(truth_f)

  de.genes <- lapply(res, function(r){
    df <- r$df
    deg <- rownames(df[which(df$padj < 0.05),])
    frac <- data.frame(
                 DE = sum(grepl("^DE", deg)),
                 DM = sum(grepl("^DM", deg)),
                 DB = sum(grepl("^DB", deg)),
                 DP = sum(grepl("^DP", deg)),
                 EE = sum(grepl("^EE", deg)),
                 EP = sum(grepl("^EP", deg)))
    frac$precision <- sum(grepl("^D", deg)) / length(deg)
    frac$recall <- sum(grepl("^D", deg)) / 2000
    frac$f1 <- -2*(frac$precision*frac$recall)/(frac$precision + frac$recall)
    return(frac)
  })
  de.df <- do.call(rbind, de.genes)
  de.df <- reshape2::melt(t(de.df))
  colnames(de.df) <- c("type", "method", "count")
  new_methods <- unlist(sapply(as.character(de.df$method),
                               function(m) {strsplit(m, "_")[[1]][2]}))
  de.df$method <- new_methods
  de.df$count <- unlist(de.df$count)
  de.df$method_type <- paste(de.df$method, de.df$type, sep = "_")
  f1 <- de.df[which(de.df$type == "f1"), ]
  f1 <- f1[order(f1$count), ]
  de.df$method <- factor(de.df$method, levels = f1$method)

  de.df <- de.df[which(de.df$type %in% c("DE", "DM", "DB", "DP", "EE", "EP")), ]
  g <- ggplot(de.df, aes(x = method,
                         y = count, fill = type,
                         group = type)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 11),
          legend.position = "top"
          ) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(limits = c(0, 750))
  ggsave(g, width = 10, height = 7, filename = "sim.svg", units = "in")
  }

plot_speed <- function(res_files, type = "sim")
{
  list.timing <- lapply(res_files, function(res_f) {
    res <- readRDS(res_f)
    timing <- t(data.frame(sapply(res, function(r)r$timing)))
    timing <- timing[, "user.self", drop = F]
    }
  )
  time.df <- do.call(rbind, list.timing)
  time.df <- data.frame(method = rownames(time.df), time = time.df[, 1])
  time.df$type <- type
  return(time.df)
  time.df <- time.df %>% dplyr::mutate(method = forcats::fct_reorder(method, time,
                                                                 fun = mean, na.rm = FALSE,
                                                                 .desc = F))
  df3 <- data_summary(time.df, varname="time",
                      groupnames=c("method"))
  df3$type <- type

}

plot_speed_graphic <- function()
{
  require(dplyr)
  require(ggplot2)
  sim.files <- list.files("/Users/bioturing/Code/SignacBenchMark/real/")
  sim.files <- sim.files[grep("^sim", sim.files)]
  d.sim <- plot_speed(file.path("/Users/bioturing/Code/SignacBenchMark/real/", sim.files), "sim")
  null.files <- list.files("/Users/bioturing/Code/SignacBenchMark/null/UMI/rds")
  d.null <- plot_speed(file.path("/Users/bioturing/Code/SignacBenchMark/null/UMI/rds", null.files), "null-UMI")
  fullength.files <- list.files("/Users/bioturing/Code/SignacBenchMark/null/full-length/rds")
  d.fullength <- plot_speed(file.path("/Users/bioturing/Code/SignacBenchMark/null/full-length/rds", fullength.files), "null-full-length")

  d <- rbind(d.sim, d.null, d.fullength)

  d$method <- unlist(sapply(as.character(d$method),
                            function(m) {strsplit(m, "_")[[1]][2]}))
  d <- d %>% dplyr::mutate(method = forcats::fct_reorder(method, time,
                                                         fun = median, na.rm = FALSE,
                                                         .desc = F))
  colors <- c('#e6194b', '#3cb44b', '#4363d8')
  g_speed <- ggplot(d, aes(x=method, y=time, color=type)) +
    geom_boxplot() +
    theme_bw() +
    scale_color_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = c(0.2, 0.8)
    ) + xlab("Method") + ylab("time(s)") + labs(color="Dataset type")
  ggsave(g_speed, filename = file.path("figures", "speed.svg"), width = 7, height = 5, units = "in", dpi = 120, device = "svg")
}

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
