args <- (commandArgs(trailingOnly = TRUE))
data_f <- args[1]
prefix <- args[2]
outdir <- args[3]
SignacBenchmark::run_FDR(data_f, prefix = prefix, out.dir = outdir)
