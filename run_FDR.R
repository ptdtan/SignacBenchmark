args <- (commandArgs(trailingOnly = TRUE))
data_f <- args[1]
prefix <- args[2]
SignacBenchmark::run_FDR(data_f, prefix = prefix)
