.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(edgeR)
library(readr)
library(magrittr)
library(tibble)

source("/hpcfs/users/a1666761/290921_trophoblast_dtu/scripts/R_code/scripts/catchSalmon_new.R")

basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu"
dir <- paste0(basedir, "/data/quant/star_salmon_quant_gffread")

samples <- list.files(path = dir)

sa_files <- file.path(dir, samples)

names(sa_files) <- samples

print("Starting import...")

SA <- catchSalmon_new(paths = sa_files,
                      verbose = TRUE)

write_rds(SA,
          paste0(basedir, "/data/import/STAR_salmon/star_salmon_object_GRCh38.rds"))
