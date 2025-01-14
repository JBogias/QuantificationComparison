
get_stager_res <- function(dtu_results) {
  
  require(magrittr)
  require(tibble)
  require(tximport)
  require(DRIMSeq)
  require(stageR)
  require(dplyr)
  require(reshape2)
  require(ggplot2)
  
  pScreen <- DRIMSeq::results(dtu_results)$pvalue %>%
    set_names(DRIMSeq::results(dtu_results)$gene_id)
  
  pScreen <- pScreen[complete.cases(pScreen)]
  
  pConfirmation <- matrix(DRIMSeq::results(
    dtu_results,
    level = "feature"
  )$pvalue, ncol = 1) %>%
    set_rownames(DRIMSeq::results(dtu_results, level = "feature")$feature_id)
  
  pConfirmation <- pConfirmation[complete.cases(pConfirmation), ] %>%
    as.matrix()
  
  tx2gene <- DRIMSeq::results(dtu_results,
                              level = "feature")[, c("feature_id", "gene_id")]
  
  stager_DTU <- stageRTx(pScreen = pScreen,
                         pConfirmation = pConfirmation,
                         pScreenAdjusted = FALSE,
                         tx2gene = tx2gene) %>%
    stageWiseAdjustment(method = "dtu",
                        alpha = 0.05,
                        allowNA = TRUE) %>%
    getAdjustedPValues(order = TRUE,
                       onlySignificantGenes = FALSE) %>%
    dplyr::rename(Gene = geneID,
                  Transcript = txID,
                  geneFDR = gene,
                  txpFDR = transcript) %>%
    as_tibble() %>%
    mutate(Transcript = str_remove(Transcript, "\\..*")) %>%
    mutate("DTU" = geneFDR < 0.05) %>%
    mutate("txpDTU" = txpFDR < 0.05)
  
  return(stager_DTU)
}

