de_glmfit <- function(y) {
  
  require(dplyr)
  require(magrittr)
  require(tibble)
  require(edgeR)
  require(limma)
  
  # Define model matrix here, group specifies the control vs knockout design
  designMatrix <- model.matrix(~ group, data = y$samples)
  
  # Change the colnames because they have "group" in front of them
  colnames(designMatrix) <- c("Intercept", "knockout")
  
  # Also need to estimate dispersions for DE testing
  y <- estimateCommonDisp(y, design = designMatrix)
  
  # Now can use glmFit
  lambda <- 1.4
  alpha <- 0.05
  
  # Perform GLM fit on the dtelist
  yTxps <- glmFit(y, design = designMatrix) %>%
    glmTreat(coef = "knockout",
             lfc = log2(lambda)) %>%
    topTags(n = Inf) %>%
    .[["table"]] %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::mutate(DE = FDR < alpha)
  
  return(yTxps)
  
}

de_lmfit <- function(y) {
  
  require(dplyr)
  require(magrittr)
  require(tibble)
  require(edgeR)
  require(limma)
  
  # Define model matrix here, group specifies the control vs knockout design
  # We don't actually need to remove the intercept because we only have two
  # groups, therefore removing the need for a contrast matrix. However, some
  # people may review my code and criticise me for not including contrasts
  # and will not accept my explanation that I simply do not need one because
  # they refuse to accept my explanations and would prefer to hold silly dogmas
  # rather than actually learn about the statistics behind it.
  designMatrix <- model.matrix(~0 + group, data = y$samples)
  
  # Change the colnames because they have "group" in front of them
  colnames(designMatrix) <- colnames(designMatrix) %>% str_remove("group")
  
  # Contrasts
  contr.matrix <- makeContrasts(
    control_vs_knockout = control-knockout, levels = colnames(designMatrix)
  )
  
  # limma voom process now
  v <- voom(y, design = designMatrix, plot = FALSE)
  vfit <- lmFit(v, design = designMatrix)
  cvfit <- contrasts.fit(vfit, contrasts = contr.matrix)
  efit <- eBayes(cvfit)
  top_tab <- topTable(efit, coef = "control_vs_knockout",
                      number = Inf, adjust = "BH") %>%
    as_tibble() %>%
    dplyr::rename("FDR" = "adj.P.Val") %>%
    dplyr::mutate(DE = FDR < 0.05 & abs(logFC) > 1)
  
  return(top_tab)
  
}