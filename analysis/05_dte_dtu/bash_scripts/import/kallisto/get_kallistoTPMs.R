.libPaths(c("/hpcfs/users/a1666761/R_Packages",.libPaths()))

getKallisto_TPMs <- function(basedir) {
  
  require(dplyr)
  require(magrittr)
  require(tibble)
  require(readr)
  
  k_dirs <- paste0(basedir, "/", list.files(basedir))
  names(k_dirs) <- list.files(basedir)
  
  for(j in 1:length(k_dirs)) {
  
    k_dirs_df <- k_dirs %>%
      as.data.frame() %>%
      rownames_to_column("names") %>%
      set_colnames(c("names", "paths"))
    
    current_k_path <- paste0(k_dirs_df$paths[j], "/abundance.tsv")
    print(j)
    print(current_k_path)
    
    k_file <- read_tsv(current_k_path)
    
    NTx <- length(k_file$tpm)
    NSamples <- length(k_dirs)
    
    k_tpms <- dplyr::select(k_file,
                            target_id,
                            tpm)
    
    if (j == 1L) {
      tpm_df <- k_tpms
    } else if (j > 1L) {
      tpm_df <- left_join(k_tpms,
                          tpm_df,
                          by = "target_id")
    }
    
  }
  tpm_df %>% set_colnames(c("transcript_id", names(k_dirs)))
}

basedir <- "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/alignment/kallisto"

kallisto_abundance <- getKallisto_TPMs(basedir)

write_csv(kallisto_abundance, "/hpcfs/users/a1666761/290921_trophoblast_dtu/data/import/kallisto/kallisto_TPMs_GRCh38.csv.gz")
