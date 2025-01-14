get_dtu_tximport <- function(tximport_counts, annotation) {
  
  require(dplyr)
  require(magrittr)
  require(tibble)
  require(tximport)
  require(DRIMSeq)
  require(ggplot2)
  
  annotation_df <- annotation %>%
    as_tibble() %>%
    dplyr::select(
      "start",
      "end",
      "length" = "width",
      "strand",
      "gene_id",
      "gene_name",
      "gene_biotype",
      "source",
      "chromosome" = "seqnames",
      "transcript_id",
      "transcript_name",
      "type"
    ) %>%
    dplyr::filter(
      !str_detect(
        chromosome,
        "^G|^H|^K"
      )
    )
  
  # Get gene-level annotations
  gene_anno <- annotation_df %>%
    dplyr::filter(
      type == "gene"
    ) %>%
    dplyr::select(
      -transcript_id,
      -transcript_name,
      -type
    ) %>%
    mutate(
      region_key = paste0(
        chromosome,
        ":",
        start,
        "-",
        end,
        ":",
        strand
      )
    )
  
  # Get transcript-level annotations
  transcript_anno <- annotation_df %>%
    dplyr::filter(
      type == "transcript"
    ) %>%
    dplyr::select(
      -gene_id,
      -type
    ) %>%
    mutate(
      region_key = paste0(
        chromosome,
        ":",
        start,
        "-",
        end,
        ":",
        strand
      )
    )
  
  tx_names <- transcript_anno %>%
    dplyr::select(transcript_id,
                  transcript_name)
  
  # Prepare counts for DTU analysis
  dtu_counts <- tximport_counts %>%
    dplyr::filter(str_detect(transcript_id, "^E")) %>%
    dplyr::mutate(transcript_id = str_remove(transcript_id, "\\..*")) %>%
    dplyr::filter(transcript_id %in% transcript_anno$transcript_id) %>%
    as.data.frame() %>% 
    inner_join(gene_txp_anno %>%
                 dplyr::select("gene_id", "transcript_id"),
               by = "transcript_id") %>%
    dplyr::rename("feature_id" = "transcript_id") %>%
    as_tibble()
  
  # Prepare the sample metadata using the names of samples - sample info was
  # found from the original publication of this data
  sirna_ko <- c("SRR13401116",
                "SRR13401117",
                "SRR13401118",
                "SRR13401119")
  
  sirna_control <- colnames(dtu_counts)[
    !colnames(dtu_counts) %in% sirna_ko
  ]
  
  ko_df <- sirna_ko %>%
    data.frame("knockout") %>%
    set_colnames(c("ID", "group"))
  
  control_df <- sirna_control %>%
    data.frame("control") %>%
    set_colnames(c("ID", "group"))
  
  public_metadata <- rbind(ko_df, control_df) %>%
    mutate(group = factor(group, levels = c("control", "knockout"))) %>%
    dplyr::filter(!ID %in% c("feature_id", "gene_id"))
  
  # Run DRIMseq method
  set.seed(123)
  
  d <- dmDSdata(
    counts = as.data.frame(dtu_counts),
    samples = as.data.frame(public_metadata) %>%
      dplyr::rename("sample_id" = "ID")
  )
  
  inds <- 1:length(d@counts)
  mean_expression <- unlist(lapply(inds, function(g){
    mean(colSums(d@counts[[g]]), na.rm = TRUE)
  }))
  
  names(mean_expression) <- names(d@counts)
  
  design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
  
  dtu_results <- dmFilter(d,
                          min_samps_gene_expr = 4,
                          min_samps_feature_expr = 4,
                          min_samps_feature_prop = 4,
                          min_feature_prop = 0.1,
                          min_gene_expr = 10,
                          min_feature_expr = 10
  ) %>%
    dmPrecision(design = design_full) %>%
    dmFit(design = design_full) %>%
    dmTest(coef = "groupknockout")
  
  return(dtu_results)
  
}
