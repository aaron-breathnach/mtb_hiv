make_figure_2 <- function(de, ids, gene_list) {
  
  tmp_1 <- de %>%
    select(-gene_set) %>%
    distinct() %>%
    inner_join(gene_list, by = "gene", relationship = "many-to-many")
  
  sig_gen <- tmp_1 %>%
    filter(q_val <= 0.05 & abs(log_fc) >= 1) %>%
    select(gene, pathway) %>%
    distinct()
  
  tmp_2 <- tmp_1 %>%
    inner_join(sig_gen, by = c("gene", "pathway")) %>%
    mutate(y = paste(pathway, "|", gene)) %>%
    inner_join(ids, by = "comparison")
  
  y_ord <- tmp_2 %>%
    filter(group == "Mtb") %>%
    arrange(pathway, log_fc)
  
  group_ord <- c("Mtb", "HIV", "Mtb/HIV vs Ctrl", "Mtb/HIV vs HIV")
  
  dat <- tmp_2 %>%
    mutate(y = factor(y, levels = y_ord$y)) %>%
    mutate(group = factor(group, levels = group_ord)) %>%
    mutate(pathway = ifelse(grepl("Complex", pathway), "Complexes I-V", pathway)) 
  
  n <- 0.05 * abs(max(dat$log_fc) - min(dat$log_fc))
  
  ann <- dat %>%
    mutate(sig = ifelse(q_val <= 0.05 & abs(log_fc) >= 1, "*", "")) %>%
    mutate(log_fc = log_fc + (sign(log_fc) * n))
  
  p <- ggplot(dat, aes(x = log_fc, y = y)) +
    facet_grid(pathway ~ group, scales = "free_y", space = "free_y") +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0)),
      labels = function(x) gsub(".* | ", "", x)
    ) +
    geom_vline(xintercept = 0) +
    geom_bar(aes(fill = log_fc), stat = "identity", colour = "black") +
    geom_text(
      data = ann,
      aes(x = log_fc, y = y, label = sig),
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
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold", angle = 0)
    ) +
    scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
    labs(
      x = "**log<sub>2</sub> fold change**",
      y = "**Gene**",
      fill = "**-log<sub>10</sub>(*p*)**"
    )
  
  ggsave("plots/figure_2.png", p, width = 10, height = 12.5, dpi = 300)
  ggsave("plots/figure_2.tiff", p, width = 10, height = 12.5, dpi = 300)
  
}
