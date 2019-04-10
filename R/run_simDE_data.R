sim_PBMC4K_1 <- function(seed, prefix = "sim_batch2", FC = 2,
                         check_points = c(1.3, 2), n.genes = 25)
{
  master_dir <- "data.1/PBMC4k"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$kmeans$clusters[3, ]

  # Randomly select two groups, groups 1 is four fold in size
  clusters1 <- which(clusters == 1)
  n <- length(clusters1)
  set.seed(seed)
  groups.1 <- sample(clusters1, round(0.2*n), replace = F)
  groups.2 <- setdiff(clusters1, groups.1)

  rowCounts <- Matrix::rowSums(mat.sp > 0)
  rowSums <- Matrix::rowSums(mat.sp)
  rowlog2Mean <- log2(rowSums/rowCounts + 1)
  # Order gene by mean
  mat.sp <- mat.sp[order(-rowlog2Mean), ]
  rowlog2Mean <- rowlog2Mean[order(-rowlog2Mean)]
  rowCounts <- Matrix::rowSums(mat.sp > 0)
  keep <- which(rowCounts > 100)
  mat.append <- as.matrix(mat.sp[which(rowCounts <= 100), ])
  mat.sp <- as.matrix(mat.sp[keep, ])
  rowlog2Mean <- rowlog2Mean[keep]

  truth <- vector("list", 2)
  # Concat the simulated matrix

  for(k in seq(1, length(check_points))){
    exp.1.low <- check_points[k]
    exp.1.high <- check_points[k] * FC

    i <- head(which(rowlog2Mean < exp.1.low), n = n.genes)
    i.low <- c(i, i - length(i))
    print(i.low)

    stopifnot(length(i.low) == 2 * n.genes && sum(i.low < 0) == 0)

    j <- tail(which(rowlog2Mean > exp.1.high), n = n.genes)
    i.high <- c(j, j + length(j))
    print(i.high)

    stopifnot(length(i.high) == 2 * n.genes && sum(i.high < 0) == 0)

    stopifnot(length(intersect(i.low, i.high)) == 0)
    truth[[k]] <- list(low = rownames(mat.sp)[i.low],
                       high = rownames(mat.sp)[i.high])
    mat.sp[truth[[k]][["low"]], groups.1] <- mat.sp[truth[[k]][["high"]], groups.1]
  }
  names(truth) <- check_points

  dir <- file.path("./", paste(prefix, "logfc", FC, "checkpoints",
                               check_points[1], check_points[2], sep = "_"))
  dir.create(dir)
  mat <- mat.sp[, c(groups.1, groups.2)]
  rownames(mat) <- rownames(mat.sp)
  real.data.1 <- seq(1, length(groups.1))
  real.data.2 <- seq(1, length(groups.2)) + length(groups.1)
  sim_PBMC4k.1 <- list(real = list(real.data.1, real.data.2),
                       null = NULL,
                       mat = mat)
  data.truth <- file.path(dir, paste0("PBMC4k_", seed, "_sim_truth.rds"))
  data.file <- file.path(dir, paste0("PBMC4k_", seed, "_sim_data.rds"))
  saveRDS(sim_PBMC4k.1, data.file)
  saveRDS(truth, data.truth)
}

sim_PBMC4K_1_subsampling <- function(rate = 0.8)
{
  master_dir <- "data.1/PBMC4k"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$kmeans$clusters[3, ]

  message("Randomly select two groups, groups 1 is four fold in size")
  clusters1 <- which(clusters == 1)
  n <- length(clusters1)
  groups.1 <- sample(clusters1, round(0.2*n), replace = F)
  groups.2 <- setdiff(clusters1, groups.1)

  genes.mean.exp <- log2(Matrix::rowMeans(mat.sp) + 1)

  # Order gene by mean
  mat.sp <- mat.sp[order(-genes.mean.exp), ]
  genes.mean.exp <- genes.mean.exp[order(-genes.mean.exp)]

  FC <- 2
  n.genes <- 25

  check_points <- c(0.5, 1.5)
  truth <- vector("list", 2)
  message("Concat the simulated matrix")

  for(k in seq(1, length(check_points))){
    exp.1.low <- check_points[k]
    exp.1.high <- check_points[k] * FC

    i <- head(which(genes.mean.exp < exp.1.low), n = n.genes)
    i.low <- c(i, i - length(i))

    j <- head(which(genes.mean.exp < exp.1.high), n = n.genes)
    i.high <- c(j, j - length(j))

    stopifnot(length(intersect(i.low, i.high)) == 0)
    truth[[k]] <- list(low = rownames(mat.sp)[i.low],
                       high = rownames(mat.sp)[i.high])
    mat.sp[truth[[k]][["low"]], groups.1] <- mat.sp[truth[[k]][["high"]], groups.1]
  }
  names(truth) <- check_points

  message("Subsampling at rate:", rate)
  mat <- as(mat.sp[, c(groups.1, groups.2)], "matrix")
  rownames(mat) <- rownames(mat.sp)
  real.data.1 <- seq(1, length(groups.1))
  real.data.2 <- seq(1, length(groups.2)) + length(groups.1)

  m <- mat[, real.data.1]
  sub_m <- vector("list", ncol(m))
  for(i in seq(1, ncol(m))){
    message(i)
    m1 <- m[, i]
    UMIs <- unlist(lapply(seq(1, length(m1)), function(i)rep(i, m1[i])))
    set.seed(1)
    retain <- sample(UMIs, rate * length(UMIs), replace = F)
    m1.new <- rep(0, length(m1))
    retain.genes <- as.numeric(names(table(retain)))
    m1.new[retain.genes] <- table(retain)
    sub_m[[i]] <- m1.new
    message(paste("Colsum old, new: ", sum(m1), sum(m1.new)))
  }
  m <- do.call(cbind, sub_m)
  mat[, real.data.1] <- m

  sim_PBMC4k.1 <- list(real = list(real.data.1, real.data.2),
                       null = NULL,
                       mat = mat)
  data.truth <- file.path(paste("data.1/PBMC4k_1_", rate, "_sim_truth.rds", ""))
  data.file <- file.path(paste("data.1/PBMC4k_1_", rate, "_sim_data.rds", ""))
  saveRDS(sim_PBMC4k.1, data.file)
  saveRDS(truth, data.truth)
}


sim_PBMC4K_1_small <- function(seed)
{
  master_dir <- "data.1/PBMC4k"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$kmeans$clusters[3, ]

  # Randomly select two groups, groups 1 is four fold in size
  clusters1 <- which(clusters == 1)
  n <- length(clusters1)

  set.seed(seed)
  n1 <- sample(seq(20, 50))
  groups.1 <- sample(clusters1, n1, replace = F)
  groups.2 <- setdiff(clusters1, groups.1)

  genes.mean.exp <- log2(Matrix::rowMeans(mat.sp) + 1)

  # Order gene by mean
  mat.sp <- mat.sp[order(-genes.mean.exp), ]
  genes.mean.exp <- genes.mean.exp[order(-genes.mean.exp)]

  FC <- 2
  n.genes <- 25

  check_points <- c(0.5, 1.5)
  truth <- vector("list", 2)
  # Concat the simulated matrix

  for(k in seq(1, length(check_points))){
    exp.1.low <- check_points[k]
    exp.1.high <- check_points[k] * FC

    i <- head(which(genes.mean.exp < exp.1.low), n = n.genes)
    i.low <- c(i, i - length(i))

    j <- head(which(genes.mean.exp < exp.1.high), n = n.genes)
    i.high <- c(j, j - length(j))

    stopifnot(length(intersect(i.low, i.high)) == 0)
    truth[[k]] <- list(low = rownames(mat.sp)[i.low],
                       high = rownames(mat.sp)[i.high])
    mat.sp[truth[[k]][["low"]], groups.1] <- mat.sp[truth[[k]][["high"]], groups.1]
  }
  names(truth) <- check_points

  mat <- as(mat.sp[, c(groups.1, groups.2)], "matrix")
  rownames(mat) <- rownames(mat.sp)
  real.data.1 <- seq(1, length(groups.1))
  real.data.2 <- seq(1, length(groups.2)) + length(groups.1)
  sim_PBMC4k.1 <- list(real = list(real.data.1, real.data.2),
                       null = NULL,
                       mat = mat)
  data.truth <- file.path(paste0("data.1/PBMC4k_small", seed, "_sim_truth.rds"))
  data.file <- file.path(paste0("data.1/PBMC4k_small", seed, "_sim_data.rds"))
  saveRDS(sim_PBMC4k.1, data.file)
  saveRDS(truth, data.truth)
  return(list(data.file, data.truth))
}
