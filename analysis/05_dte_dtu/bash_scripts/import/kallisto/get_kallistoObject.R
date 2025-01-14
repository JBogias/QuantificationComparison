.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(edgeR)
library(readr)

basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu"
dir <- paste0(basedir, "/data/alignment/kallisto_global")

samples <- list.files(path = dir)

dirs <- file.path(dir, samples)

print("Starting import...")

kallisto <- catchKallisto(paths = dirs,
                      verbose = TRUE)

write_rds(kallisto,
          paste0(basedir, "/data/import/kallisto_object_HTR8_grch38_global.rds"))
