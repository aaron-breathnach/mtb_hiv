make_heatmap <- function(pwy, genes, norm_data, metadata, gois) {
  
  ann_col <- metadata %>%
    select(sample_id, group) %>%
    rename(Group = 2) %>%
    column_to_rownames("sample_id")
  
  ann_row <- gois %>%
    select(1, 3) %>%
    rename(Complex = 2) %>%
    arrange(Complex) %>%
    column_to_rownames("gene")
  
  complexes <- ann_row %>%
    pull(Complex) %>%
    unique()
  
  pal <- c(
    "Ctrl" = "#bbbbbb",
    "TB" = "#0096ff",
    "TB/HIV" = "#9437ff",
    "HIV" = "#ff2600"
  )
  
  col_complexes <- c("#D7191C", "#FDAE61", "#FFFFBF", "#ABDDA4", "#2B83BA")
  names(col_complexes) <- complexes
  
  ann_colours <- list(
    "Group" = pal,
    "Complex" = col_complexes
  )
  
  mat <- norm_data %>%
    filter(gene %in% genes) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  mat <- mat[rownames(ann_row),]
  
  log_mat <- log(mat + 1)
  
  newnames <- lapply(
    rownames(log_mat),
    function(x) bquote(italic(.(x))))
  
  p <- pheatmap::pheatmap(
    log_mat,
    scale = "row",
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_colors = ann_colours,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    labels_row = as.expression(newnames)
  )
  
  ggplotify::as.ggplot(p)
  
}

make_figure_4  <- function(norm_data, metadata, gene_list) {
  
  gois <- gene_list %>%
    filter(grepl("Complex", pathway)) %>%
    arrange(subpathway, gene)
  
  genes <- gois$gene
  
  p_a <- viz_pathway("Complexes I-V", gois, de, "heatmap")
  p_b <- make_heatmap(pwy, genes, norm_data, metadata, gois)
  p_c <- viz_pathway("Complexes I-V", gois, de, "barchart")
  
  top <- cowplot::plot_grid(
    plotlist = list(plot.new(), p_a, plot.new()),
    rel_widths = c(1, 2, 1),
    nrow = 1
  )
  
  left <- cowplot::plot_grid(
    plotlist = list(top, p_c),
    rel_heights = c(1, 0.75),
    nrow = 2,
    scale = 0.9,
    labels = c("A", "C")
  )
  
  right <- cowplot::plot_grid(
    plotlist = list(plot.new(), p_b, plot.new()),
    rel_heights = c(1, 8, 1),
    nrow = 3,
    labels = c("", "B")
  )
  
  p <- cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(2, 1))
  
  ggsave("plots/figure_4.png", p, width = 20, height = 15, dpi = 300, bg = "white")
  ggsave("plots/figure_4.tiff", p, width = 20, height = 15, dpi = 300, bg = "white")
  
}
