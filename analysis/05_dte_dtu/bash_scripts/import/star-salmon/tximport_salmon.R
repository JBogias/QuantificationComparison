.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(tximport)
library(readr)
library(magrittr)
library(tibble)


basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/quant/star_salmon_quant"

sa_files <- paste0(basedir, "/", list.files(basedir), "/quant.sf")

names(sa_files) <- list.files(basedir)

transCounts_sa <- tximport(sa_files,
                           type = "salmon",
                           txOut = TRUE,
                           countsFromAbundance = "scaledTPM")

cts <- transCounts_sa$counts[rowSums(transCounts_sa$counts) > 0,]

cts %>%
   as.data.frame() %>%
   rownames_to_column("transcript_id") %>%
   as_tibble() %>%
   write_csv("/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/STAR_salmon/cts_star_salmon.csv.gz")

write_rds(transCounts_sa,
          "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/STAR_salmon/tximport_star_salmon.rds")
