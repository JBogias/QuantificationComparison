
plot_ma <- function(x) {
  
  require(dplyr)
  require(ggplot2)
  require(magrittr)
  require(tibble)
  require(edgeR)
  require(limma)
  
  x %>%
  mutate(group = case_when(
    FDR < 0.05 & logFC > 1 ~ "Up Regulated",
    FDR < 0.05 & logFC < -1 ~ "Down Regulated",
    .default = "Not Significant"
  ),
  group = factor(group, levels = c("Up Regulated", "Down Regulated",
                                   "Not Significant"))) %>%
  ggplot(aes(x = AveExpr, y = logFC, colour = group)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("red", "royalblue", "grey80")) +
  labs(x = "Mean Transcript Expression (log2 CPM)",
       y = "Log2 Fold Change",
       colour = "") +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 13))
}

plot_volc <- function(x) {
  
  require(dplyr)
  require(ggplot2)
  require(magrittr)
  require(tibble)
  require(edgeR)
  require(limma)
  
  x %>%
    mutate(group = case_when(
      FDR < 0.05 & logFC > 1 ~ "Up Regulated",
      FDR < 0.05 & logFC < -1 ~ "Down Regulated",
      .default = "Not Significant"
    ),
    group = factor(group, levels = c("Up Regulated", "Down Regulated",
                                     "Not Significant"))) %>%
    ggplot(aes(x = logFC, y = -log10(FDR), colour = group)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
    geom_vline(xintercept = -1, linetype = "dashed", colour = "blue") +
    scale_color_manual(values = c("red", "royalblue", "grey80")) +
    labs(x = "Log2 Fold Change",
         y = "False Discovery Rate (-log10)",
         colour = "") +
    theme_bw() +
    theme(axis.title = element_text(colour = "black", size = 14),
          axis.text = element_text(colour = "black", size = 12),
          legend.text = element_text(colour = "black", size = 13))
}