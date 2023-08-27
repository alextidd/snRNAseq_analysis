# config
base_dir <- ifelse(Sys.info()["login"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "~/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

# load packages
library(SCPA)
library(msigdbr)
library(magrittr)
library(dplyr)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(circlize)
library(tidyverse)

# functions
fix_pathway_names <- function(pathway_names) {
  pathway_names %>% tolower() %>% gsub("\\_", " ", .) %>% gsub("hallmark ", "", .)
}
plot_my_heatmap <- function(scpa_out_df, col, stat = "qval") {
  scpa_out_df %>%
    select(Pathway, name = matches(col), value = matches(stat)) %>%
    mutate(Pathway = fix_pathway_names(Pathway)) %>%
    pivot_wider() %>%
    column_to_rownames("Pathway") %>%
    as.matrix() %>%
    Heatmap(name = stat,
            border = T,
            show_row_dend = F,
            show_column_dend = T,
            row_names_gp = grid::gpar(fontsize = 8))
}
plot_my_enrichment <- function(p_dat_in, top_n = 5) {
  p_dat <-
    p_dat_in %>%
    mutate(Pathway = Pathway %>% fix_pathway_names(),
           color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                             FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                             FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                             FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

  ggplot(p_dat, aes(-FC, qval)) +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
    geom_point(cex = 2.6, shape = 21, fill = p_dat$color, stroke = 0.3) +
    ggrepel::geom_text_repel(data = p_dat %>% slice_max(abs(FC), n = top_n), aes(label = Pathway),
                             size = 7, min.segment.length = 0, nudge_x = 10) +
    xlab("Enrichment") +
    ylab("Qval") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          aspect.ratio = 1) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
}

# load cds
seu <-
  readRDS("output/snRNAseq_workflow/by_patient_wo_organoids/O/integrating/infercnv/malignant_integrating/seu.rds")

# get hallmark pathways
pathways <-
  msigdbr(species = "Homo sapiens", category = "H") %>%
  format_pathways()



