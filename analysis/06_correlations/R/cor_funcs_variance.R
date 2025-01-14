get_joined_samples <- function(w, x, y, z,
                               method_w = "method_w", method_x = "method_x",
                               method_y = "method_y", method_z = "method_z",
                               sample = 1) {
  
  require(dplyr)
  require(magrittr)
  require(tibble)
  require(ggplot2)
  require(edgeR)
  
  cpm_to_df <- function(x) {
    xdf <- x %>%
      cpm(log = TRUE) %>%
      as.data.frame() %>%
      mutate(mean_exp = rowMeans(cpm(x, log = TRUE))) %>%
      rownames_to_column("transcript_id") 
    
    return(xdf)
  }
  
  sample_no <- colnames(x)[sample]
  
  wdf <- cpm_to_df(w)
  colnames(wdf) <- c("transcript_id", paste0(method_w, "_", colnames(wdf)[-1]))
  
  xdf <- cpm_to_df(x)
  colnames(xdf) <- c("transcript_id", paste0(method_x, "_", colnames(xdf)[-1]))
  
  ydf <- cpm_to_df(y)
  colnames(ydf) <- c("transcript_id", paste0(method_y, "_", colnames(ydf)[-1]))
  
  zdf <- cpm_to_df(z)
  colnames(zdf) <- c("transcript_id", paste0(method_z, "_", colnames(zdf)[-1]))
  
  jdf <- inner_join(wdf, xdf,
                    by = "transcript_id") %>%
    inner_join(ydf,
               by = "transcript_id") %>%
    inner_join(zdf,
               by = "transcript_id") %>%
    as_tibble()
  
  jdf %>%
    dplyr::select(
      "transcript_id",
      paste0(method_w, "_", sample_no),
      paste0(method_x, "_", sample_no),
      paste0(method_y, "_", sample_no),
      paste0(method_z, "_", sample_no)
    )
}

# Plot sample correlations ----
plot_sample_cor <- function(df, x, y, variance_df) {
  
  require(dplyr)
  require(magrittr)
  require(tibble)
  require(ggplot2)
  
  cohort_x <- colnames(df)[x]
  cohort_y <- colnames(df)[y]
  
  method_x <- colnames(df)[x] %>% str_remove("_SRR.*")
  method_y <- colnames(df)[y] %>% str_remove("_SRR.*")
  
  sample_id <- colnames(df)[x] %>% str_remove(".*_")
  
  concordant <- head(dplyr::arrange(variance_df, 
                                    rowvar)[["transcript_id"]], 100)
  discordant <- head(dplyr::arrange(variance_df,
                                    desc(rowvar))[["transcript_id"]], 100)
  
  cor_to_plot <- df %>%
    #dplyr::arrange(colour_txps) %>%
    dplyr::mutate(colour_txps = case_when(
      transcript_id %in% concordant ~ "Concordant",
      transcript_id %in% discordant ~ "Discordant",
      .default = "midcor"
    ),
    colour = factor(colour_txps, levels = c("midcor",
                                            "Concordant",
                                            "Discordant")
    ))
  
  suppressWarnings(cor_legend <- cowplot::get_legend(
    cor_to_plot %>%
      dplyr::filter(!colour_txps == "midcor") %>%
      ggplot() +
      geom_point(
        aes(x = .data[[cohort_x]],
            y = .data[[cohort_y]],
            colour = .data[["colour_txps"]]),
        alpha = 1,
        size = 8
      ) +
      scale_colour_manual(
        values = c("green", "red", "black")
      ) +
      guides(colour = guide_legend(bycol = TRUE)) +
      labs(colour = "Transcript Counts\nVariance") +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.text = element_text(size = 20, colour = "black"),
            legend.spacing.x = unit(1, "cm"),
            legend.position = "bottom")
  ))
  
  
  cor_colours <- cor_to_plot %>%
    dplyr::filter(colour_txps %in% c("Concordant", "Discordant"))
  col_scale <- c("green", "red")
  
  
  cor_to_plot %>%
    ggplot() +
    geom_point(
      aes(x = .data[[cohort_x]],
          y = .data[[cohort_y]]),
      colour = "black",
      alpha = 1
    ) +
    geom_point(data = cor_colours,
               aes(
                 x = .data[[cohort_x]],
                 y = .data[[cohort_y]],
                 colour = .data[["colour_txps"]])
    ) +
    scale_colour_manual(values = col_scale) +
    geom_smooth(
      aes(x = .data[[cohort_x]], y = .data[[cohort_y]]),
      method = "gam",
      linewidth = 0.5
    ) +
    labs(x = paste0(method_x),
         y = paste0(method_y)) +
    theme_bw() +
    theme(axis.title = element_text(colour = "black",
                                    size = 18),
          axis.text = element_text(colour = "black",
                                   size = 14),
          legend.position = "none")
}

# Get Correlation legend ----
get_cor_legend <- function(df, x, y, variance_df) {
  require(tidyverse)
  require(ggplot2)
  
  cohort_x <- colnames(df)[x]
  cohort_y <- colnames(df)[y]
  
  method_x <- colnames(df)[x] %>% str_remove("_SRR.*")
  method_y <- colnames(df)[y] %>% str_remove("_SRR.*")
  
  sample_id <- colnames(df)[x] %>% str_remove(".*_")
  
  concordant <- head(dplyr::arrange(variance_df,
                                    rowvar)[["transcript_id"]], 100)
  discordant <- head(dplyr::arrange(variance_df,
                                    desc(rowvar))[["transcript_id"]], 100)
  
  cor_to_plot <- df %>%
    #dplyr::arrange(colour_txps) %>%
    dplyr::mutate(colour_txps = case_when(
      transcript_id %in% concordant ~ "Concordant",
      transcript_id %in% discordant ~ "Discordant",
      .default = "midcor"
    ),
    colour = factor(colour_txps, levels = c("midcor",
                                            "Concordant",
                                            "Discordant")
    ))
  
  cor_legend <- cowplot::get_legend(
    cor_to_plot %>%
      dplyr::filter(colour_txps %in% c("Concordant", "Discordant")) %>%
      ggplot() +
      geom_point(
        aes(x = .data[[cohort_x]],
            y = .data[[cohort_y]],
            colour = .data[["colour_txps"]]),
        alpha = 1,
        size = 8
      ) +
      scale_colour_manual(
        values = c("green", "red")
      ) +
      guides(colour = guide_legend(bycol = TRUE)) +
      labs(colour = "Variance") +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.text = element_text(size = 20, colour = "black"),
            legend.spacing.x = unit(1, "cm"),
            legend.position = "bottom")
  )
  return(cor_legend)
}

# Plot the correlation values ----
plot_cor_value <- function(r_sq) {
  
  require(tidyverse)
  require(ggplot2)
  require(magrittr)
  require(edgeR)
  
  ggplot() + 
    annotate("text",
             x = 1, 
             y = 1, 
             size = 20,
             label = round(r_sq, 3)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black",
                                      size = 1,
                                      fill = NA),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

# Plot the entire correlation cowplot ----
plot_cor_cowplot_var <- function(joined_df, variance_df) {
  
  require(tidyverse)
  require(ggplot2)
  require(grid)
  require(magrittr)
  require(cowplot)
  require(edgeR)
  
  # Now plot all and then create cowplot
  a_b <- plot_sample_cor(df = joined_df, x = 2, y = 3,
                         variance_df = variance_df) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  a_b_cor <- cor(joined_df[[2]],
                 joined_df[[3]]) %>%
    plot_cor_value() +
    labs(x = "Kallisto", y = "") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  
  a_c <- plot_sample_cor(df = joined_df, x = 2, y = 4,
                         variance_df = variance_df) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  a_c_cor <- cor(joined_df[[2]],
                 joined_df[[4]]) %>%
    plot_cor_value() +
    labs(x = "STAR", y = "") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  a_d <- plot_sample_cor(df = joined_df, x = 2, y = 5,
                         variance_df = variance_df)
  a_d_cor <- cor(joined_df[[2]],
                 joined_df[[5]]) %>%
    plot_cor_value() +
    labs(x = "Salmon", y = "HISAT2") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  b_c <- plot_sample_cor(df = joined_df, x = 3, y = 4,
                         variance_df = variance_df) +
    theme(axis.title = element_blank(),
          axis.text = element_blank())
  
  b_c_cor <- cor(joined_df[[3]],
                 joined_df[[4]]) %>%
    plot_cor_value() +
    labs(x = "", y = "")
  
  b_d <- plot_sample_cor(df = joined_df, x = 3, y = 5,
                         variance_df = variance_df) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  b_d_cor <- cor(joined_df[[3]],
                 joined_df[[5]]) %>%
    plot_cor_value() +
    labs(x = "", y = "Kallisto") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  c_d <- plot_sample_cor(df = joined_df, x = 4, y = 5,
                         variance_df = variance_df) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  c_d_cor <- cor(joined_df[[4]],
                 joined_df[[5]]) %>%
    plot_cor_value() +
    labs(x = "", y = "STAR") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  blanc <- ggplot() +
    geom_blank() +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  cor_cowplot <- cowplot::plot_grid(blanc, a_b_cor, a_c_cor, a_d_cor,
                                    a_b, blanc, b_c_cor, b_d_cor,
                                    a_c, b_c, blanc, c_d_cor,
                                    a_d, b_d, c_d, blanc,
                                    ncol = 4,
                                    nrow = 4,
                                    rel_heights = c(1,1,1,1),
                                    rel_widths = c(1,1,1,1))
  
  xlab <- grid::textGrob("Transcript Expression (log2 CPM)",
                         gp = gpar(fontface = "bold",
                                   col = "black",
                                   fontsize = 26))
  
  fig_cor <- cowplot::plot_grid(cor_cowplot, xlab,
                                ncol = 1, rel_heights = c(1, 0.05))
  
  ylab <- grid::textGrob("Transcript Expression (log2 CPM)",
                         gp = gpar(fontface = "bold",
                                   col = "black",
                                   fontsize = 26),
                         rot = 90)
  
  cor_legend <- get_cor_legend(df = joined_df, x = 2, y = 3, 
                               variance_df = variance_df)
  
  fig_cor_labs <- cowplot::plot_grid(ylab, fig_cor,
                                     blanc, cor_legend,
                                     ncol = 2,
                                     nrow = 2,
                                     rel_heights = c(1, 0.05),
                                     rel_widths = c(0.05, 1))
  return(fig_cor_labs)
}

plot_cor_cowplot_loadings <- function(joined_df, colour_txps = NA) {
  
  require(dplyr)
  require(tibble)
  require(ggplot2)
  require(grid)
  require(magrittr)
  require(cowplot)
  require(edgeR)
  
  # First define sample plotting function
  plot_sample_cor <- function(df, x, y) {
    
    require(tidyverse)
    require(ggplot2)
    
    cohort_x <- colnames(df)[x]
    cohort_y <- colnames(df)[y]
    
    method_x <- colnames(df)[x] %>% str_remove("_SRR.*")
    method_y <- colnames(df)[y] %>% str_remove("_SRR.*")
    
    sample_id <- colnames(df)[x] %>% str_remove(".*_")
    
    cor_to_plot <- df %>%
      #dplyr::arrange(colour_txps) %>%
      dplyr::mutate(colour_txps = case_when(
        transcript_id %in% head(hisat2_loadings$transcript_id, 100) ~ "Top 100 Negative\nPC1 Loadings",
        transcript_id %in% head(star_loadings$transcript_id, 100) ~ "Top 100 Negative\nPC2 Loadings",
        transcript_id %in% head(salmon_loadings$transcript_id, 100) ~ "Top 100 Positive\nPC2 Loadings",
        transcript_id %in% head(kallisto_loadings$transcript_id, 100) ~ "Top 100 Positive\nPC6 Loadings",
        .default = "Not a Loading"
      ),
      colour = factor(colour_txps, levels = c("Not a Loading",
                                              "Top 100 Negative\nPC1 Loadings",
                                              "Top 100 Negative\nPC2 Loadings",
                                              "Top 100 Positive\nPC2 Loadings",
                                              "Top 100 Positive\nPC6 Loadings")
      ))
    
    cor_legend <- cowplot::get_legend(
      cor_to_plot %>%
        ggplot() +
        geom_point(
          aes(x = .data[[cohort_x]],
              y = .data[[cohort_y]],
              colour = .data[["colour_txps"]]),
          alpha = 1,
          size = 8
        ) +
        scale_colour_manual(
          values = c("#000000", "#61c3d7", "#4fc14d", "#bf4dc1", "#f6866f")
        ) +
        guides(colour = guide_legend(bycol = TRUE)) +
        labs(colour = "Loadings") +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 20, colour = "black"),
              legend.spacing.x = unit(1, "cm"),
              legend.position = "bottom")
    )
    
    if (x == 2 & y == 3) {
      cor_colours <- cor_to_plot %>%
        dplyr::filter(!colour_txps %in% c("Not a Loading",
                                          "Top 100 Positive\nPC2 Loadings",
                                          "Top 100 Positive\nPC6 Loadings"))
      col_scale <- c("#61c3d7", "#4fc14d")
    } else if (x == 2 & y == 4) {
      cor_colours <- cor_to_plot %>%
        dplyr::filter(!colour_txps %in% c("Not a Loading",
                                          "Top 100 Negative\nPC2 Loadings",
                                          "Top 100 Positive\nPC6 Loadings"))
      col_scale <- c("#61c3d7", "#bf4dc1")
    } else if (x == 2 & y == 5) {
      cor_colours <- cor_to_plot %>%
        dplyr::filter(!colour_txps %in% c("Not a Loading",
                                          "Top 100 Negative\nPC2 Loadings",
                                          "Top 100 Positive\nPC2 Loadings"))
      col_scale <- c("#61c3d7", "#f6866f")
    } else if (x == 3 & y == 4) {
      cor_colours <- cor_to_plot %>%
        dplyr::filter(!colour_txps %in% c("Not a Loading",
                                          "Top 100 Negative\nPC1 Loadings",
                                          "Top 100 Positive\nPC6 Loadings"))
      col_scale <- c("#4fc14d", "#bf4dc1")
    } else if (x == 3 & y == 5) {
      cor_colours <- cor_to_plot %>%
        dplyr::filter(!colour_txps %in% c("Not a Loading",
                                          "Top 100 Negative\nPC1 Loadings",
                                          "Top 100 Positive\nPC2 Loadings"))
      col_scale <- c("#4fc14d", "#f6866f")
    } else if (x == 4 & y == 5) {
      cor_colours <- cor_to_plot %>%
        dplyr::filter(!colour_txps %in% c("Not a Loading",
                                          "Top 100 Negative\nPC1 Loadings",
                                          "Top 100 Negative\nPC2 Loadings"))
      col_scale <- c("#bf4dc1", "#f6866f")
    }
    
    cor_to_plot %>%
      ggplot() +
      geom_point(
        aes(x = .data[[cohort_x]],
            y = .data[[cohort_y]]),
        colour = "black",
        alpha = 1
      ) +
      geom_point(data = cor_colours,
                 aes(
                   x = .data[[cohort_x]],
                   y = .data[[cohort_y]],
                   colour = .data[["colour_txps"]])
      ) +
      #scale_colour_manual(values = c("grey80","#dedb0d", "#43bf46", "#4e63df", "#f6866f")) +
      scale_colour_manual(values = col_scale) +
      geom_smooth(
        aes(x = .data[[cohort_x]], y = .data[[cohort_y]]),
        method = "gam"
      ) +
      labs(x = paste0(method_x),
           y = paste0(method_y)) +
      theme_bw() +
      theme(axis.title = element_text(colour = "black",
                                      size = 18),
            axis.text = element_text(colour = "black",
                                     size = 14),
            legend.position = "none")
  }
  
  get_cor_legend <- function(df, x, y) {
    require(tidyverse)
    require(ggplot2)
    
    cohort_x <- colnames(df)[x]
    cohort_y <- colnames(df)[y]
    
    method_x <- colnames(df)[x] %>% str_remove("_SRR.*")
    method_y <- colnames(df)[y] %>% str_remove("_SRR.*")
    
    sample_id <- colnames(df)[x] %>% str_remove(".*_")
    
    cor_to_plot <- df %>%
      #dplyr::arrange(colour_txps) %>%
      dplyr::mutate(colour_txps = case_when(
        transcript_id %in% head(hisat2_loadings$transcript_id, 100) ~ "Top 100 Negative\nPC1 Loadings",
        transcript_id %in% head(kallisto_loadings$transcript_id, 100) ~ "Top 100 Negative\nPC2 Loadings",
        transcript_id %in% head(star_loadings$transcript_id, 100) ~ "Top 100 Positive\nPC2 Loadings",
        transcript_id %in% head(salmon_loadings$transcript_id, 100) ~ "Top 100 Positive\nPC6 Loadings",
        .default = "Not a Loading"
      ),
      colour = factor(colour_txps, levels = c("Not a Loading",
                                              "Top 100 Negative\nPC1 Loadings",
                                              "Top 100 Negative\nPC2 Loadings",
                                              "Top 100 Positive\nPC2 Loadings",
                                              "Top 100 Positive\nPC6 Loadings")
      ))
    
    cor_legend <- cowplot::get_legend(
      cor_to_plot %>%
        ggplot() +
        geom_point(
          aes(x = .data[[cohort_x]],
              y = .data[[cohort_y]],
              colour = .data[["colour_txps"]]),
          alpha = 1,
          size = 8
        ) +
        scale_colour_manual(
          values = c("#000000", "#61c3d7", "#4fc14d", "#bf4dc1", "#f6866f")
        ) +
        guides(colour = guide_legend(bycol = TRUE)) +
        labs(colour = "Loadings") +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 20, colour = "black"),
              legend.spacing.x = unit(1, "cm"),
              legend.position = "bottom")
    )
    return(cor_legend)
  }
  
  # Now define correlation annotation function
  plot_cor_value <- function(r_sq) {
    
    require(tidyverse)
    require(ggplot2)
    require(magrittr)
    require(edgeR)
    
    ggplot() + 
      annotate("text",
               x = 1, 
               y = 1, 
               size = 20,
               label = round(r_sq, 3)) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(position = "right") +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black",
                                        size = 1,
                                        fill = NA),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  }
  
  # if (is.na(colour_txps[1])) {
  #   joined_df$colour <- "Yes"
  # } else {
  #   joined_df <- joined_df %>%
  #     mutate(
  #       colour = if_else(transcript_id %in% colour_txps,
  #                        "Yes",
  #                        "No"
  #       )
  #     )
  # }
  
  # Now plot all and then create cowplot
  a_b <- plot_sample_cor(df = joined_df, x = 2, y = 3) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  a_b_cor <- cor(joined_df[[2]],
                 joined_df[[3]]) %>%
    plot_cor_value() +
    labs(x = "Kallisto", y = "") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  
  a_c <- plot_sample_cor(df = joined_df, x = 2, y = 4) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  a_c_cor <- cor(joined_df[[2]],
                 joined_df[[4]]) %>%
    plot_cor_value() +
    labs(x = "STAR", y = "") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  a_d <- plot_sample_cor(df = joined_df, x = 2, y = 5)
  a_d_cor <- cor(joined_df[[2]],
                 joined_df[[5]]) %>%
    plot_cor_value() +
    labs(x = "Salmon", y = "HISAT2") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  b_c <- plot_sample_cor(df = joined_df, x = 3, y = 4) +
    theme(axis.title = element_blank(),
          axis.text = element_blank())
  
  b_c_cor <- cor(joined_df[[3]],
                 joined_df[[4]]) %>%
    plot_cor_value() +
    labs(x = "", y = "")
  
  b_d <- plot_sample_cor(df = joined_df, x = 3, y = 5) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  b_d_cor <- cor(joined_df[[3]],
                 joined_df[[5]]) %>%
    plot_cor_value() +
    labs(x = "", y = "Kallisto") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  c_d <- plot_sample_cor(df = joined_df, x = 4, y = 5) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  c_d_cor <- cor(joined_df[[4]],
                 joined_df[[5]]) %>%
    plot_cor_value() +
    labs(x = "", y = "STAR") +
    theme(axis.title = element_text(colour = "black",
                                    size = 18))
  
  blanc <- ggplot() +
    geom_blank() +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  cor_cowplot <- cowplot::plot_grid(blanc, a_b_cor, a_c_cor, a_d_cor,
                                    a_b, blanc, b_c_cor, b_d_cor,
                                    a_c, b_c, blanc, c_d_cor,
                                    a_d, b_d, c_d, blanc,
                                    ncol = 4,
                                    nrow = 4,
                                    rel_heights = c(1,1,1,1),
                                    rel_widths = c(1,1,1,1))
  
  xlab <- grid::textGrob("Transcript Expression (log2 CPM)",
                         gp = gpar(fontface = "bold",
                                   col = "black",
                                   fontsize = 26))
  
  fig_cor <- cowplot::plot_grid(cor_cowplot, xlab,
                                ncol = 1, rel_heights = c(1, 0.05))
  
  ylab <- grid::textGrob("Transcript Expression (log2 CPM)",
                         gp = gpar(fontface = "bold",
                                   col = "black",
                                   fontsize = 26),
                         rot = 90)
  
  cor_legend <- get_cor_legend(df = joined_df, x = 1, y = 2)
  
  fig_cor_labs <- cowplot::plot_grid(ylab, fig_cor,
                                     blanc, cor_legend,
                                     ncol = 2,
                                     nrow = 2,
                                     rel_heights = c(1, 0.05),
                                     rel_widths = c(0.05, 1))
}
  
  # Get the facet correlations ----
  get_facet_cor <- function(x, y, 
                            method_x = "method_x",
                            method_y = "method_y") {
    
    require(tidyverse)
    require(ggplot2)
    require(magrittr)
    require(edgeR)
    
    cpm_to_df <- function(x) {
      xdf <- x %>%
        cpm(log = TRUE) %>%
        as.data.frame() %>%
        mutate(`Mean Expression` = rowMeans(cpm(x, log = TRUE))) %>%
        rownames_to_column("transcript_id") 
      
      return(xdf)
    }
    
    xdf <- cpm_to_df(x)
    colnames(xdf) <- c("transcript_id", paste0(method_x, ":", colnames(xdf)[-1]))
    
    ydf <- cpm_to_df(y)
    colnames(ydf) <- c("transcript_id", paste0(method_y, ":", colnames(ydf)[-1]))
    
    jdf <- inner_join(xdf, ydf,
                      by = "transcript_id") %>%
      pivot_longer(cols = starts_with(c("b", "s", "k")),
                   names_to = "id",
                   values_to = "expr") %>%
      separate(col = id, into = c("aligner", "sample"), sep = ":") %>%
      as_tibble()
    
    jdfx <- jdf %>%
      dplyr::filter(aligner == method_x) %>%
      dplyr::rename("aligner_x" = "aligner",
                    "expr_x" = "expr")
    
    jdfy <- jdf %>%
      dplyr::filter(aligner == method_y) %>%
      dplyr::rename("aligner_y" = "aligner",
                    "expr_y" = "expr")
    
    inner_join(jdfx, jdfy,
               by = c("transcript_id", "sample")) %>%
      ggplot() +
      geom_point(aes(x = expr_x, y = expr_y), alpha = 0.2, size = 1) +
      geom_smooth(aes(x = expr_x, y = expr_y), method = "lm") +
      labs(x = paste0(method_x, " Transcript Expression (log2 CPM)"),
           y = paste0(method_y, " Transcript Expression (log2 CPM)")) +
      facet_wrap(~ sample) +
      theme_bw()
  }
  
  
  