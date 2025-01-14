
draw_venn <- function(hisat2_res,
                      kallisto_res,
                      star_res,
                      salmon_res) {
  
  # Load required packages
  require(dplyr)
  require(tidyr)
  require(magrittr)
  require(tibble)
  require(VennDiagram)
  
  # Set names for significance column to know if certain transcripts were
  # DE or not
  
  if("DE" %in% colnames(hisat2_res)) {
    
    hisat2_glmfit_fdrs <- hisat2_res$DE %>%
      set_names(hisat2_res$transcript_id)
    
    kallisto_glmfit_fdrs <- kallisto_res$DE %>%
      set_names(kallisto_res$transcript_id)
    
    star_glmfit_fdrs <- star_res$DE %>%
      set_names(star_res$transcript_id)
    
    salmon_glmfit_fdrs <- salmon_res$DE %>%
      set_names(salmon_res$transcript_id)
    
  } else if("DTU" %in% colnames(hisat2_res)) {
    
    hisat2_glmfit_fdrs <- hisat2_res$DTU %>%
      set_names(hisat2_res$Transcript)
    
    kallisto_glmfit_fdrs <- kallisto_res$DTU %>%
      set_names(kallisto_res$Transcript)
    
    star_glmfit_fdrs <- star_res$DTU %>%
      set_names(star_res$Transcript)
    
    salmon_glmfit_fdrs <- salmon_res$DTU %>%
      set_names(salmon_res$Transcript)
    
  } else {
    
    hisat2_glmfit_fdrs <- hisat2_res$txpDTU %>%
      set_names(hisat2_res$Transcript)
    
    kallisto_glmfit_fdrs <- kallisto_res$txpDTU %>%
      set_names(kallisto_res$Transcript)
    
    star_glmfit_fdrs <- star_res$txpDTU %>%
      set_names(star_res$Transcript)
    
    salmon_glmfit_fdrs <- salmon_res$txpDTU %>%
      set_names(salmon_res$Transcript)
    
  }
  
  # Define the areas for the venn diagram
  area1 <- sum(hisat2_glmfit_fdrs)
  area2 <- sum(kallisto_glmfit_fdrs)
  area3 <- sum(star_glmfit_fdrs)
  area4 <- sum(salmon_glmfit_fdrs)
  
  # Need to extract the names for intersection to plot venn diagram
  hisat2_txps <- names(hisat2_glmfit_fdrs[hisat2_glmfit_fdrs])
  kallisto_txps <- names(kallisto_glmfit_fdrs[kallisto_glmfit_fdrs])
  star_txps <- names(star_glmfit_fdrs[star_glmfit_fdrs])
  salmon_txps <- names(salmon_glmfit_fdrs[salmon_glmfit_fdrs])
  
  # Define intersects to be used as plotting data
  n12 <- length(intersect(hisat2_txps, kallisto_txps))
  n13 <- length(intersect(hisat2_txps, star_txps))
  n14 <- length(intersect(hisat2_txps, salmon_txps))
  
  n23 <- length(intersect(kallisto_txps, star_txps))
  n24 <- length(intersect(kallisto_txps, salmon_txps))
  
  n34 <- length(intersect(star_txps, salmon_txps))
  
  n123 <- length(
    Reduce(intersect, list(hisat2_txps, kallisto_txps, star_txps))
  )
  n124 <- length(
    Reduce(intersect, list(hisat2_txps, kallisto_txps, salmon_txps))
  )
  
  n134 <- length(
    Reduce(intersect, list(hisat2_txps, star_txps, salmon_txps))
  )
  
  n234 <- length(
    Reduce(intersect, list(kallisto_txps, star_txps, salmon_txps))
  )
  
  n1234 <- length(
    Reduce(intersect, list(hisat2_txps, kallisto_txps, star_txps, salmon_txps))
  )
  
  # Plot out venn diagram with defined objects
  venn.plot <- draw.quad.venn(
    area1 = area1, area2 = area2, area3 = area3, area4 = area4,
    n12 = n12, n13 = n13, n14 = n14,
    n23 = n23, n24 = n24,
    n34 = n34,
    n123 = n123, n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = c("HISAT2", "Kallisto", "STAR", "Salmon"),
    fill = c("red", "blue", "yellow", "green"),
    cat.col = c("darkred", "navy", "darkorange1", "darkgreen"),
    ind = TRUE)
}

