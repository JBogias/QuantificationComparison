.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(edgeR)
library(readr)

source("/hpcfs/users/a1666761/290921_trophoblast_dtu/scripts/R_code/catchSalmon_new.R")

basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu"
dir <- paste0(basedir, "/data/alignment/selective_alignment")

samples <- list.files(path = dir)

dirs <- file.path(dir, samples)

print("Starting import...")

SA <- catchSalmon_new(paths = dirs,
                      verbose = TRUE)

write_rds(SA,
          paste0(basedir, "/data/import/SA-SalmonObject_HTR8_grch38.rds"))
