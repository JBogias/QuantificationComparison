filterCounts <- function (x, minCPM = 2, minSamp = 4) {
  require(edgeR)
  
  cpm <- edgeR::cpm(x, log = FALSE)
  i <- rowSums(cpm > minCPM) > minSamp
  x[i,]
}

get_filter_plot <- function(dat_object,
                            bootstrapped = "none",
                            xlab_off = TRUE,
                            labels = c("a", "b"),
                            min_count = 1,
                            min_sample = 4) {
  
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(tximport)
  require(ggplot2)
  require(magrittr)
  require(tibble)
  require(edgeR)
  require(limma)
  
  # Creating a DGEList object for use in edgeR.
  # The lightweight mappers kallisto, salmon, and SA can all utilise the 
  # overdispersion estimates generated from the alignment stage of kallisto and
  # salmon
  
  if (tolower(bootstrapped) == "bootstrapped") {
    
    dte_counts <- dat_object %>%
      with(counts / annotation$Overdispersion) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      dplyr::filter(str_detect(transcript_id, "^E")) %>%
      as.data.frame() %>%
      column_to_rownames("transcript_id") %>%
      as.matrix() %>%
      set_colnames(basename(colnames(.))) %>%
      set_rownames(str_remove(rownames(.), "\\..*")) %>%
      cpm(log = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      pivot_longer(cols = starts_with("SRR"),
                   names_to = "sample_id",
                   values_to = "logcpm")
    
    dte_counts_filtered <- dat_object %>%
      with(counts / annotation$Overdispersion) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      dplyr::filter(str_detect(transcript_id, "^E")) %>%
      as.data.frame() %>%
      column_to_rownames("transcript_id") %>%
      as.matrix() %>%
      filterCounts(minCPM = min_count, minSamp = min_sample) %>%
      set_colnames(basename(colnames(.))) %>%
      set_rownames(str_remove(rownames(.), "\\..*")) %>%
      cpm(log = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      pivot_longer(cols = starts_with("SRR"),
                   names_to = "sample_id",
                   values_to = "logcpm")
    
  } else if (tolower(bootstrapped) == "none") {
    
    dte_counts <- dat_object$counts %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      dplyr::filter(str_detect(transcript_id, "^E")) %>%
      as.data.frame() %>%
      column_to_rownames("transcript_id") %>%
      as.matrix() %>%
      set_colnames(basename(colnames(.))) %>%
      set_rownames(str_remove(rownames(.), "\\..*")) %>%
      cpm(log = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      pivot_longer(cols = starts_with("SRR"),
                   names_to = "sample_id",
                   values_to = "logcpm")
    
    dte_counts_filtered <- dat_object$counts %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      dplyr::filter(str_detect(transcript_id, "^E")) %>%
      as.data.frame() %>%
      column_to_rownames("transcript_id") %>%
      as.matrix() %>%
      filterCounts(minCPM = min_count, minSamp = min_sample) %>%
      set_colnames(basename(colnames(.))) %>%
      set_rownames(str_remove(rownames(.), "\\..*")) %>%
      cpm(log = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      pivot_longer(cols = starts_with("SRR"),
                   names_to = "sample_id",
                   values_to = "logcpm")
    
  } else {
    message(
      paste0("Not valid method. Either pick bootstrapped or none")
    )
  }
  
  unfiltered_plot <- dte_counts %>%
    ggplot() +
    geom_density(aes(x = logcpm, colour = sample_id)) +
    labs(x = "Transcript Expression (log2 CPM)",
         y = "Density",
         colour = "Sample ID") +
    { if(xlab_off) labs(x = "") } +
    coord_cartesian(xlim = c(-4, 12.5),
                    ylim = c(0, 2)) +
    theme_bw() +
    theme(axis.title = element_text(colour = "black",
                                    size = 12),
          axis.text = element_text(colour = "black",
                                   size = 10),
          legend.position = "none")
  
  filtered_plot <- dte_counts_filtered %>%
    ggplot() +
    geom_density(aes(x = logcpm, colour = sample_id)) +
    labs(x = "Transcript Expression (log2 CPM)",
         y = "",
         colour = "Sample ID") +
    { if(xlab_off) labs(x = "") } +
    coord_cartesian(xlim = c(-4, 12.5),
                    ylim = c(0, 0.5)) +
    theme_bw() +
    theme(axis.title = element_text(colour = "black",
                                    size = 12),
          axis.text = element_text(colour = "black",
                                   size = 10),
          legend.title = element_text(colour = "black",
                                      size = 11),
          legend.text = element_text(colour = "black",
                                     size = 8))
  
  filplot_legend <- cowplot::get_legend(filtered_plot)
  
  filplot_legend %>%
    write_rds(here("data/cowplot_filter_legend.rds"))
  
  filtered_plot <- filtered_plot +
    theme(legend.position = "none")
  
  filter_cowplot <- cowplot::plot_grid(unfiltered_plot,
                                       filtered_plot,
                                       nrow = 1,
                                       labels = labels,
                                       rel_widths = c(1, 1))
  
  filter_cowplot
  
}
