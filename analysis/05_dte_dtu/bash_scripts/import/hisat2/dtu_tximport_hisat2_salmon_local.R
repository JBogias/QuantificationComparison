.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(tximport)
library(readr)
library(magrittr)
library(tibble)
library(here)

#tx2gene <- read_csv("/hpcfs/users/a1666761/Refs/ref_annotations/tx2gene_grch38_103.csv.gz")
tx2gene <- read_csv(here("data/annotations/txp_gene_ensdb_lengths.csv.gz")) %>%
  dplyr::select("tx_id_version", "gene_id")

#basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/quant/hisat2_salmon"
basedir <- here("data/counts/hisat2_quant")

hisat2_files <- paste0(basedir, "/", list.files(basedir), "/quant.sf")

names(hisat2_files) <- list.files(basedir)

transCounts_hisat2 <- tximport(hisat2_files,
                               type = "salmon",
                               txOut = TRUE,
                               tx2gene = tx2gene,
                               countsFromAbundance = "dtuScaledTPM")

cts <- transCounts_hisat2$counts[rowSums(transCounts_hisat2$counts) > 0,]

cts %>%
   as.data.frame() %>%
   rownames_to_column("transcript_id") %>%
   as_tibble() %>%
   write_csv(here("data/counts/dtu_cts_hisat2_salmon.csv.gz"))

write_rds(transCounts_hisat2,
          "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/hisat2_salmon/dtu_tximport_hisat2_salmon.rds")
