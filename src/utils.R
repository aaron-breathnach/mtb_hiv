get_dat <- function(path) {
  
  comparison <- path %>%
    basename() %>%
    str_replace(".de_res.tsv", "")
  
  cols <- c("Log2 fold change", "P-value", "BH.p.value", "Gene.sets")
  
  read_delim(path) %>%
    drop_na() %>%
    select(1, all_of(cols)) %>%
    setNames(c("gene", "log_fc", "p_val", "q_val", "gene_set")) %>%
    mutate(gene = gsub("-mRNA", "", gene)) %>%
    separate_rows(gene_set, sep = ", ") %>%
    mutate(comparison = comparison)
  
}

viz_pathway <- function(pwy, gois, de, method) {
  
  ids <- get_ids()
  
  gois <- gois %>%
    filter(pathway == pwy)
  
  dat <- de %>%
    select(-gene_set) %>%
    distinct() %>%
    inner_join(gois, by = "gene") %>%
    inner_join(ids, by = "comparison")
  
  n <- 0.05 * abs(max(dat$log_fc) - min(dat$log_fc))
  
  n_gene <- dat %>%
    pull(gene) %>%
    unique() %>%
    length()
  
  if (pwy == "Complexes I-V") {
    
    p_inp <- dat %>%
      select(gene, subpathway, log_fc, group) %>%
      pivot_wider(names_from = group, values_from = log_fc, values_fill = 0) %>%
      pivot_longer(!c(gene, subpathway), names_to = "group", values_to = "log_fc")
    
  } else {
    
    p_inp <- dat %>%
      select(gene, log_fc, group) %>%
      pivot_wider(names_from = group, values_from = log_fc, values_fill = 0) %>%
      pivot_longer(!gene, names_to = "group", values_to = "log_fc")
    
  }
  
  group_ord <- c("Mtb", "HIV", "Mtb/HIV vs Ctrl", "Mtb/HIV vs HIV")
  
  p_inp$group <- factor(p_inp$group, levels = group_ord)
  
  gene_ord <- p_inp %>%
    filter(group == "Mtb") %>%
    arrange(log_fc) %>%
    pull(gene)
  
  p_inp$gene <- factor(p_inp$gene, levels = gene_ord)
  
  ann <- dat %>%
    mutate(sig = ifelse(q_val <= 0.05, "*", "")) %>%
    mutate(log_fc = log_fc + (sign(log_fc) * n)) %>%
    mutate(group = factor(group, levels = group_ord)) %>%
    mutate(gene = factor(gene, levels = gene_ord))
  
  if (method == "heatmap") {
    
    p <- ggplot(p_inp, aes(x = group, y = gene)) +
      geom_tile(aes(fill = log_fc), colour = "black") +
      geom_text(
        data = ann,
        aes(label = sig),
        vjust = 0.75,
        size = 5,
        colour = "black"
      ) +
      scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      scale_x_discrete(expand = expansion(mult = c(0, 0))) +
      scale_y_discrete(expand = expansion(mult = c(0, 0)), position = "right") +
      labs(x = "Group", y = "Gene", fill = "**log<sub>2</sub> fold change**") +
      theme_bw(base_size = 12.5) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(face = "bold"),
        legend.title = ggtext::element_markdown(),
        strip.text.y.left = element_text(face = "bold", angle = 0)
      )
    
    if (pwy == "Complexes I-V") {
      p <- p + 
        facet_grid(
          subpathway ~ .,
          scales = "free",
          space = "free",
          switch = "both"
        )
    }
    
  } else {
    
    p <- ggplot(p_inp, aes(x = log_fc, y = gene)) +
      facet_grid(~ group, scales = "free_y", space = "free_y") +
      scale_y_discrete(
        expand = expansion(mult = c(0, 0)),
        labels = function(x) gsub(".* | ", "", x)
      ) +
      geom_vline(xintercept = 0) +
      geom_bar(aes(fill = log_fc), stat = "identity", colour = "black") +
      geom_text(
        data = ann,
        aes(x = log_fc, y = gene, label = sig),
        vjust = 0.75,
        size = 7.5
      ) +
      theme_bw(base_size = 12.5) +
      theme(
        axis.text.y = element_text(face = "italic"),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown(),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold", angle = 0)
      ) +
      scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      labs(
        x = "**log<sub>2</sub> fold change**",
        y = "**Gene**",
        fill = "**log<sub>2</sub> fold change**"
      )
  }
  
  return(p)
  
}

get_ids <- function() {
  tibble(comparison = c(
    "U1_Infected_vs_U1_Uninfected",
    "U1_Infected_vs_U937_Uninfected",
    "U1_Uninfected_vs_U937_Uninfected",
    "U937_Infected_vs_U937_Uninfected"
  ),
  group = c(
    "Mtb/HIV vs HIV",
    "Mtb/HIV vs Ctrl",
    "HIV",
    "Mtb"
  ))
}
