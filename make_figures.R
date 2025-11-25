library(tidyverse)

source("src/make_figure_1.R")
source("src/make_figure_2.R")
source("src/make_figure_3.R")
source("src/make_figure_4.R")
source("src/utils.R")

make_figures <- function() {
  
  if (!dir.exists("plots")) dir.create("plots")
  
  norm_data <- read_delim("data/norm_data.tsv")
  
  metadata <- read_tsv("data/metadata.tsv") %>%
    mutate(sample_id = as.character(sample_id)) %>%
    dplyr::rename(comparison = group) %>%
    mutate(group = case_when(
      hiv == "negative" & tb == "negative" ~ "Ctrl",
      hiv == "negative" & tb == "positive" ~ "TB",
      hiv == "positive" & tb == "negative" ~ "HIV",
      hiv == "positive" & tb == "positive" ~ "TB/HIV"
    )) %>%
    arrange(group)
  
  de <- list.files("data/nsolver", recursive = TRUE, full.names = TRUE) %>%
    purrr::map(\(x) get_dat(x)) %>%
    bind_rows()
  
  gene_list <- read_delim("data/gois.txt") %>%
    mutate(subpathway = pathway) %>%
    mutate(pathway = ifelse(grepl("Complex", pathway), "Complexes I-V", pathway))
  
  make_figure_1(norm_data, metadata, de)
  make_figure_2(de, ids = get_ids(), gene_list)
  make_figure_3(norm_data, metadata, gene_list)
  make_figure_4(norm_data, metadata, gene_list)
  
}
