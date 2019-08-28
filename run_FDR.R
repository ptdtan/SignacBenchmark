args <- (commandArgs(trailingOnly = TRUE))
data_f <- args[1]
outdir <- args[2]
filtering <- args[3]
SignacBenchmark::run_FDR(data_f, prefix = basename(data_f), out.dir = outdir,
                         filtering = filtering == "TRUE")
