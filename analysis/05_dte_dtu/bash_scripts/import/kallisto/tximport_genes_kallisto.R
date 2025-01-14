.libPaths(c("/hpcfs/users/a1666761/R_Packages", .libPaths()))

library(tximport)
library(readr)
library(magrittr)
library(tibble)

tx2gene <- read_csv("/hpcfs/users/a1666761/Refs/ref_annotations/tx2gene_grch38_103.csv.gz")

basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/alignment/kallisto"

sa_files <- paste0(basedir, "/", list.files(basedir), "/quant.sf")

names(sa_files) <- list.files(basedir)

transCounts_sa <- tximport(sa_files,
                               type = "salmon",
                               txOut = FALSE,
			       tx2gene = tx2gene,
                               countsFromAbundance = "scaledTPM")

cts <- transCounts_sa$counts[rowSums(transCounts_sa$counts) > 0,]

cts %>%
   as.data.frame() %>%
   rownames_to_column("transcript_id") %>%
   as_tibble() %>%
   write_csv("/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/cts_genes_kallisto.csv.gz")

write_rds(transCounts_sa,
          "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/tximport_genes_kallisto.rds")
