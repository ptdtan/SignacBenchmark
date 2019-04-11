#!/bin/bash

LIST_DATASET_FILES=$1
RSCRIPT="R_MAX_NUM_DLLS=1000 R_LIBS=/raidixshare01/single_cell_data/gene_marker/pbmc4k/Seurat/R-3.5.3/R/lib/R/library /raidixshare01/single_cell_data/gene_marker/pbmc4k/Seurat/R-3.5.3/bin/R"
cat $LIST_DATASET_FILES | parallel -j 2 "${RSCRIPT} ./run_Islam.R data.real/{}.data.rds data.real/{}.truth.rds {}"
