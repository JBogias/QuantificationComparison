.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(tximport)
library(tidyverse)
library(readr)
library(magrittr)
library(tibble)


basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/quant/star_salmon_quant"
tx2gene <- read_csv("/hpcfs/users/a1666761/Refs/ref_annotations/tx2gene_grch38_103.csv.gz")

tx2gene <- dplyr::mutate(tx2gene,
			 tx_id = str_remove(tx2gene$tx_id_version, "\\..*"),
			 tx_id_version = NULL)

tx2gene <- dplyr::relocate(tx2gene, tx_id, gene_id)

sa_files <- paste0(basedir, "/", list.files(basedir), "/quant.sf")

names(sa_files) <- list.files(basedir)

transCounts_sa <- tximport(sa_files,
                           type = "salmon",
                           txOut = TRUE,
			   tx2gene = tx2gene,
                           countsFromAbundance = "dtuScaledTPM")

cts <- transCounts_sa$counts[rowSums(transCounts_sa$counts) > 0,]

cts %>%
   as.data.frame() %>%
   rownames_to_column("transcript_id") %>%
   as_tibble() %>%
   write_csv("/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/STAR_salmon/dtu_cts_star_salmon.csv.gz")

write_rds(transCounts_sa,
          "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/STAR_salmon/dtu_tximport_star_salmon.rds")
