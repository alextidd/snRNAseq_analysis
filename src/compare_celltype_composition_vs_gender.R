base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "~/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

library(tidyverse)
library(magrittr)
library(ggpubr)
library(rstatix)
library(purrr)

# load
singler_annots <-
  readr::read_tsv("output/celltype_composition/infercnv_singler_annots.tsv")

# get celltype composition within samples
props <-
  singler_annots %>%
  #filter(lineage == "immune") %>%
  dplyr::select(patient, gender, sample_site, timepoint, label.main, label.fine, label.infercnv) %>%
  mutate(label.fine = case_when(label.infercnv == "Malignant_cells" ~ "Malignant_cells",
                                .default = label.fine),
         label.main = case_when(label.infercnv == "Malignant_cells" ~ "Malignant_cells",
                                .default = label.main),
         label.tcells = case_when(grepl("T_cell:CD8+", label.infercnv) ~ "T_cells:CD8+",
                                  grepl("T_cell:CD4+", label.infercnv) ~ "T_cells:CD4+",
                                  .default = NA)) %>%
  pivot_longer(cols = c("label.main", "label.fine", "label.infercnv")) %>%
  group_by(patient, gender, sample_site, timepoint, name, value) %>%
  summarise(n = n()) %>%
  group_by(patient, gender, sample_site, timepoint, name) %>%
  mutate(prop = n / sum(n)) %>%
  group_by(name, value, gender)

# check summary stats
props %>%
  get_summary_stats(prop, type = "mean_sd")
  
# run t-test
props %>%
  group_by(name, value, gender) %>%
  filter(n() > 1) %>%
  group_by(name, value) %>%
  filter(all(c("F", "M") %in% gender)) %>%
  t_test(prop ~ gender) %>%
  add_significance() %>% 
  arrange(p)

