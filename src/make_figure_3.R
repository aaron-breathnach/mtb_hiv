make_figure_3  <- function(norm_data, metadata, gene_list) {
  
  gois <- gene_list %>%
    filter(pathway == "Glycolysis")
  
  genes <- gois$gene
  
  p_a <- viz_pathway("Glycolysis", gois, de, "heatmap")
  p_b <- viz_pathway("Glycolysis", gois, de, "barchart")
  
  p <- cowplot::plot_grid(
    plotlist = list(p_a, p_b),
    nrow = 1,
    rel_widths = c(1, 2),
    axis = "tb",
    align = "h",
    labels = "AUTO"
  )
  
  ggsave(
    "plots/figure_3.png",
    p,
    width = 12.5,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    "plots/figure_3.tiff",
    p,
    width = 12.5,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
}
