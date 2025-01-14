theme_justin <- function() {
  require(ggplot2)
  theme(axis.title = element_text(colour = "black",
                                  size = 12),
        axis.text = element_text(colour = "black",
                                 size = 10),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey95"))
}

theme_justin_facets <- function() {
  require(ggplot2)
  theme(axis.title = element_text(colour = "black",
                                  size = 12),
        axis.text = element_text(colour = "black",
                                 size = 10),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.1),
        strip.background = element_rect(fill = "lightblue",
                                        colour = "black"),
        plot.background = element_rect(colour = "black",
                                       linewidth = 0.5),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.5),
        panel.border = element_rect(colour = "black",
                                    linewidth = 0.5,
                                    fill = NA),
        panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey95"))
}