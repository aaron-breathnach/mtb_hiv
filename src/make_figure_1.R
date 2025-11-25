#########################
## Figure 1A: PCA plot ##
#########################

make_figure_1a <- function(norm_data, metadata, pal) {
  
  pca_inp <- norm_data %>%
    pivot_longer(!gene) %>%
    pivot_wider(names_from = gene, values_from = value) %>%
    column_to_rownames("name")
  
  pca <- prcomp(pca_inp, scale. = TRUE)
  
  pca_summary <- summary(pca)$importance[2,]
  pc1_label <- paste0("PC1 (", 100 * round(pca_summary[[1]], 2), "%)")
  pc2_label <- paste0("PC2 (", 100 * round(pca_summary[[2]], 2), "%)")
  
  pca_points <- pca$x[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    mutate(sample_id = as.character(sample_id)) %>%
    inner_join(metadata, by = "sample_id")
  
  ggplot(pca_points, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = group), pch = 21, size = 3) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(face = "bold")) +
    scale_fill_manual(values = pal) +
    labs(x = pc1_label, y = pc2_label, fill = "Group")
  
}

###################################
## Figure 1B: heatmap of all DEG ##
###################################

make_figure_1b <- function(norm_data, metadata, pal, genes) {
  
  annotation <- metadata %>%
    select(sample_id, group) %>%
    rename(Group = 2) %>%
    column_to_rownames("sample_id")
  
  ann_colours <- list("Group" = pal)
  
  norm_data_sig <- norm_data %>%
    filter(gene %in% genes)
  
  mat <- norm_data_sig %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  newnames <- lapply(
    rownames(mat),
    function(x) bquote(italic(.(x)))
  )
  
  log_mat <- log(mat + 1)
  
  p <- pheatmap::pheatmap(
    log_mat,
    scale = "row",
    annotation_col = annotation,
    annotation_colors = ann_colours,
    show_colnames = FALSE,
    labels_row = as.expression(newnames)
  )
  
  ggplotify::as.ggplot(p)
  
}

##############################
## Figure 1C: volcano plots ##
##############################

volcano_plot <- function(df, title = NULL, q = 0.05) {
  
  tab <- df %>%
    mutate(colour = case_when(
      padj >  0.05 & abs(log2FoldChange) < 1 ~ "NS",
      padj >  0.05 & abs(log2FoldChange) > 1 ~ "log<sub>2</sub> FC",
      padj <= 0.05 & abs(log2FoldChange) < 1 ~ "*p*-value",
      padj <= 0.05 & abs(log2FoldChange) > 1 ~ "*p*-value and log<sub>2</sub> FC"
    ))
  
  x_min <- min(tab$log2FoldChange)
  x_max <- max(tab$log2FoldChange)
  x_val <- max(abs(c(x_min, x_max)))
  x_lim <- c(-x_val, x_val)
  
  ann <- tab %>%
    drop_na() %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
    group_by(direction) %>%
    top_n(-5, padj)
  
  pal <- c("grey", "#ff5100", "#f00080", "#7d5df9")
  
  names(pal) <- c(
    "NS",
    "log<sub>2</sub> FC",
    "*p*-value",
    "*p*-value and log<sub>2</sub> FC"
  )
  
  tmp_max_sig <- tab %>%
    filter(padj <= q)
  
  tmp_volcano <- ggplot(tab, aes(x = log2FoldChange, y = -log10(pvalue)))
  
  if (nrow(tmp_max_sig) > q) {
    
    max_sig <- filter(tmp_max_sig, pvalue == max(pvalue)) %>%
      .[[1, "pvalue"]]
    
    tmp_volcano <- tmp_volcano +
      geom_hline(yintercept = -log10(max_sig), linetype = "dashed", colour = "darkgrey")
    
  }
  
  volcano <- tmp_volcano +
    geom_vline(xintercept = c(-2, -1, 1, 2), linetype = "dashed", colour = "darkgrey") +
    geom_point(aes(colour = colour, fill = colour), pch = 21, alpha = 0.75) +
    theme_bw(base_size = 12.5) +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = ggtext::element_markdown()) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    ggtitle(title) +
    labs(x = "**log<sub>2</sub> fold change**",
         y = "**-log<sub>10</sub>(*p*)**") +
    xlim(x_lim) +
    ggrepel::geom_text_repel(
      data = ann,
      aes(x = log2FoldChange, y = -log10(pvalue), label = gene)
    )
  
  return(volcano)
  
}

wrapper_volcano <- function(x, title = NULL) {
  
  print(x)
  
  comparison <- x %>%
    str_split("/") %>%
    unlist() %>%
    nth(3)
  
  if (comparison == "U937_Infected_vs_U937_Uninfected") {
    title <- "Mtb"
  } else if (comparison == "U1_Uninfected_vs_U937_Uninfected") {
    title <- "HIV"
  } else {
    title <- "Mtb/HIV"
  }
  
  cols <- c("Log2 fold change", "P-value", "BH.p.value")
  
  tab <- read_delim(x) %>%
    select(1, all_of(cols)) %>%
    setNames(c("gene", "log2FoldChange", "pvalue", "padj")) %>%
    mutate(gene = gsub("-.*", "", gene))
  
  p <- volcano_plot(tab, title = title, q = 0.05)
  
  return(p)
  
}

make_figure_1c <- function() {
  
  filenames <- sprintf(
    "data/nsolver/%s.de_res.tsv",
    c(
      "U937_Infected_vs_U937_Uninfected",
      "U1_Uninfected_vs_U937_Uninfected",
      "U1_Infected_vs_U937_Uninfected"
    )
  )
  
  plot_list <- lapply(filenames, wrapper_volcano)
  
  plot_list[[2]] <- plot_list[[2]] +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(override.aes = list(size = 3))) 
  
  cowplot::plot_grid(
    plotlist = plot_list,
    nrow = 1,
    scale = 0.95,
    align = "h",
    axis = "tb"
  )
  
}

##############
## Figure 1 ##
##############

make_figure_1 <- function(norm_data, metadata, de) {
  
  pal <- c(
    "Ctrl" = "#bbbbbb",
    "TB" = "#0096ff",
    "TB/HIV" = "#9437ff",
    "HIV" = "#ff2600"
  )
  
  genes <- de %>%
    filter(q_val < 0.05 & abs(log_fc) > 2) %>%
    pull(gene) %>%
    unique()
  
  p_a <- make_figure_1a(norm_data, metadata, pal)
  p_b <- make_figure_1b(norm_data, metadata, pal, genes)
  p_c <- make_figure_1c()
  
  par(mar = c(1, 1, 1, 1))
  
  p <- cowplot::plot_grid(
    plotlist = list(plot.new(), p_a, plot.new()),
    nrow = 1,
    rel_widths = c(1, 1, 1),
    labels = c("", "A", "")
  )
  
  p <- cowplot::plot_grid(
    plotlist = list(p, p_c),
    nrow = 2,
    align = "v",
    axis = "lr",
    labels = c("", "C")
  )
  
  p <- cowplot::plot_grid(
    plotlist = list(p, p_b),
    nrow = 1,
    rel_widths = c(1, 1/3),
    labels = c("", "B")
  )
  
  ggsave(
    "plots/figure_1.png",
    p,
    width = 20,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    "plots/figure_1.tiff",
    p,
    width = 20,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
}
