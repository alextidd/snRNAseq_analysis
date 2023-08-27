base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "~/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

library(tidyverse)
library(magrittr)
library(ggplot2)
library(mlbench)
library(MASS)
library(pROC)
library(GGally)

# load sample metadata
sample_metadata <-
  readr::read_tsv("../cellranger/output/sample_metadata.tsv") %>%
  mutate(treatment = case_when(treatment == "NONE" ~ "unknown", TRUE ~ treatment),
         broad_response = case_when(response %in% c("Responder", "Partial Responder") ~ "responder",
                                    response %in% c("Non-responder", "Progressor") ~ "nonresponder") %>%
           factor(levels = c("nonresponder", "responder")))

# get patient-level variables
patient_metadata <-
  sample_metadata %>%
  distinct(patient_id, broad_response, gender, age, treatment) %>%
  mutate(across(where(is.character), as.factor))

# load label.infercnv
annots <-
  readr::read_tsv("output/celltype_composition/infercnv_singler_annots.tsv")

# get celltype props
celltype_props <- 
  annots %>%
  dplyr::count(patient_id, label.infercnv) %>%
  group_by(patient_id) %>% 
  transmute(label.infercnv, prop = n / sum(n)) %>% 
  pivot_wider(
    names_from = "label.infercnv", values_from = "prop", 
    names_prefix = "prop_celltype_", values_fill = 0)

# get lineage props
lineage_props <-
  annots %>%
  mutate(
    lineage = case_when(
      label.infercnv == "Malignant_cells" ~ "malignant",
      label.infercnv %in% c("NK_cells", "Myeloid_cells", "Erythrocytes", 
                            "Macrophage", "T_cells", "B_cells") ~ "immune",
      .default = "stromal")) %>%
  dplyr::count(patient_id, lineage) %>%
  group_by(patient_id) %>% 
  transmute(lineage, prop = n / sum(n)) %>% 
  pivot_wider(
    names_from = "lineage", values_from = "prop", 
    names_prefix = "prop_lineage_", values_fill = 0)

# malignancy scores in malignant cells
tumour_malignancy <-
  annots %>%
  filter(label.infercnv == "Malignant_cells") %>%
  group_by(patient_id) %>%
  summarise(tumour_malignancy = mean(malignancy))

# cycling scores in malignant cells
pse <- readRDS("output/glmGamPoi/pse.rds")
tumour_cycling <-
  colData(pse) %>%
  as_tibble() %>%
  filter(label.infercnv == "Malignant_cells") %>%
  group_by(patient_id) %>%
  summarise(tumour_cycling = mean(cycling_score))

# get final model data
model_data <-
  patient_metadata %>% 
  inner_join(celltype_props) %>%
  inner_join(lineage_props) %>%
  inner_join(tumour_malignancy) %>%
  inner_join(tumour_cycling)
summary(model_data)

# run logistic regression 
# -> binary dependent variable = broad_response
model <-
  lme4::glmer(
    broad_response ~ gender + age + 
      prop_lineage_malignant + prop_lineage_immune + prop_lineage_stromal +
      tumour_malignancy + tumour_cycling +
      (1 | treatment),
    data = model_data,
    family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa"),
    nAGQ = 10
  )

print(model, corr = F)


