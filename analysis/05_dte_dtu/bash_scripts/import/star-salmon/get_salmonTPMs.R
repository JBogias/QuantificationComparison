.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(edgeR)
library(readr)
library(magrittr)
library(tibble)

source("/hpcfs/users/a1666761/210316_GDM_RNAseq/scripts/R/catchSalmon_TPMs.R")

basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu"
dir <- paste0(basedir, "/data/quant/star_salmon_quant_gffread")

samples <- list.files(path = dir)

dirs <- file.path(dir, samples)

print("Starting import...")

SA <- catchSalmon_new(paths = dirs,
                      verbose = TRUE)

write_rds(SA,
          paste0(basedir, "/data/import/STAR_salmon/star_salmon_abundance_GRCh38.rds"))
