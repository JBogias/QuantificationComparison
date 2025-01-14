filterCounts <- function (x, minCPM = 2, minSamp = 4) {
  require(edgeR)
  cpm <- edgeR::cpm(x, log = FALSE)
  i <- rowSums(cpm > minCPM) > minSamp
  x[i,]
}

get_dtelist <- function(dat_object, bootstrapped = FALSE,
                        min_counts = 2, min_samples = 4,
                        annotations = "auto",
                        normalise = TRUE) {
  require(dplyr)
  require(tidyr)
  require(readr)
  require(tximport)
  require(magrittr)
  require(tibble)
  require(edgeR)
  require(limma)
  # Creating a DGEList object for use in edgeR.
  # The lightweight mappers kallisto, salmon, and SA can all utilise the 
  # overdispersion estimates generated from the alignment stage of kallisto and
  # salmon
  if (bootstrapped == TRUE) {
    message("Starting Bootstrapping process...")
    dte_counts <- dat_object %>%
      with(counts / annotation$Overdispersion) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      dplyr::filter(str_detect(transcript_id, "^E")) %>%
      as.data.frame() %>%
      column_to_rownames("transcript_id") %>%
      as.matrix() %>%
      filterCounts(minCPM = min_counts, minSamp = min_samples) %>%
      set_colnames(basename(colnames(.))) %>%
      set_rownames(str_remove(rownames(.), "\\..*"))
  } else if (bootstrapped == FALSE) {
    message("Bootstrapping deactivated, running without it")
    dte_counts <- dat_object$counts %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      dplyr::filter(str_detect(transcript_id, "^E")) %>%
      as.data.frame() %>%
      column_to_rownames("transcript_id") %>%
      as.matrix() %>%
      filterCounts(minCPM = min_counts, minSamp = min_samples) %>%
      set_colnames(basename(colnames(.))) %>%
      set_rownames(str_remove(rownames(.), "\\..*"))
  } else {
    message(
      paste0("Is it bootstrapped? TRUE or FALSE? Binary decision here.")
    )
  }
  # Create design matrix by using the names of samples
  # I went by what was described in the paper itself where this data is from
  sirna_ko <- c("SRR13401116",
                "SRR13401117",
                "SRR13401118",
                "SRR13401119")
  
  sirna_control <- colnames(dte_counts)[
    !colnames(dte_counts) %in% sirna_ko
  ]
  ko_df <- sirna_ko %>%
    data.frame("knockout") %>%
    set_colnames(c("ID", "group"))
  control_df <- sirna_control %>%
    data.frame("control") %>%
    set_colnames(c("ID", "group"))
  public_metadata <- rbind(ko_df, control_df) %>%
    mutate(group = factor(group, levels = c("control", "knockout")))
  # Can create the DGElist now. Simply input the counts we made before and
  # Use TMM normalisation. I don't think tximport can be normalised by TMMs
  if(isTRUE(identical(annotations, "auto")) && normalise == TRUE) {
    message(paste0("Automatic annotations specified,",
                   " attempting to load data from data/ repo"))
    tx_key <- read_csv(here("data/grch38_103_tx_name_key.csv.gz"))
    tx_anno <- read_csv(
      here("data/txp_gene_ensdb_lengths.csv.gz")
    ) %>%
      dplyr::mutate(transcript_id = tx_id) %>%
      as.data.frame() %>%
      left_join(tx_key, by = "transcript_id") %>%
      column_to_rownames("tx_id")
    y <- DGEList(dte_counts,
                 samples = public_metadata,
                 genes = tx_anno[rownames(dte_counts), ]) %>%
      calcNormFactors(method = "TMM")
  } else if(isTRUE(identical(annotations, "auto")) && normalise == FALSE) {
    message(paste0("Automatic annotations specified,",
                   " attempting to load data from data/ repo.",
                   "No normalisation set, outputting raw counts"))
    tx_key <- read_csv(here("data/grch38_103_tx_name_key.csv.gz"))
    tx_anno <- read_csv(
      here("data/txp_gene_ensdb_lengths.csv.gz")
    ) %>%
      dplyr::mutate(transcript_id = tx_id) %>%
      as.data.frame() %>%
      column_to_rownames("tx_id") %>%
      left_join(tx_key, by = "transcript_id") # nolint
    y <- DGEList(dte_counts,
                 samples = public_metadata,
                 genes = tx_anno[rownames(dte_counts), ])
  } else if(isFALSE(identical(annotations, "auto")) && normalise == TRUE) {
    message("Making DGEList with preset annotations and normalisating")
    tx_anno <- dplyr::filter(annotations,
                             transcript_id %in% rownames(dte_counts))
    
    y <- DGEList(dte_counts,
                 samples = public_metadata,
                 genes = tx_anno) %>%
      calcNormFactors(method = "TMM")
  } else {
    message("Making DGEList with preset annotations and no normalisation")
    tx_anno <- dplyr::filter(annotations,
                             transcript_id %in% rownames(dte_counts))
    
    y <- DGEList(dte_counts,
                 samples = public_metadata,
                 genes = tx_anno)
  }
  return(y)
}
