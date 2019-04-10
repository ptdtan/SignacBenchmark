args <- (commandArgs(trailingOnly = TRUE))

data_f <- args[1]
data_truth <- args[2]
prefix <- args[3]

HuyMarkerBenchmark::run_TPR_real(data_f, data_truth, prefix = prefix)
