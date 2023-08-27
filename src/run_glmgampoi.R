base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "~/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

library(glmGamPoi)
library(tidyverse)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(SEtools)
library(ComplexHeatmap)

# function: volcano plot
plot_volcano <- function(p_dat, title = "", subtitle = "") {
  dat <- p_dat %>%
    mutate(perturbation = case_when(adj_pval > 0.01 ~ NA,
                                    lfc > 1.5 ~ 'upregulated',
                                    lfc < -1.5 ~ 'downregulated')) %>%
    group_by(perturbation) %>%
    mutate(n_degs = ifelse(is.na(perturbation), NA, paste(n(), perturbation))) 
  dat %>%
    ggplot(aes(x = lfc, y = -log10(adj_pval), label = name, colour = perturbation)) +
    geom_vline(xintercept = c(-1.5, 0, 1.5), linetype = "dashed", colour = "grey") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey") +
    geom_point(size = 0.6, alpha = 0.8) +
    geom_text_repel(
      data = function(p_dat) {
        p_dat %>% 
          group_by(perturbation) %>% 
          filter(abs(lfc) > 1.5, adj_pval < 0.01) %>%
          mutate(rank_lfc = rank(-abs(lfc)), rank_adj_pval = rank(adj_pval)) %>%
          filter(rank_lfc <= 3 | rank_adj_pval <= 5)
      }, min.segment.length = 0) +
    theme_classic() +
    theme(legend.position = "none", 
          text = element_text(size = 12),
          axis.line = element_blank(), axis.ticks = element_blank()) +
    scale_colour_manual(values = c("blue", "red")) +
    scale_x_continuous(breaks = c(-1.5, 0, 1.5)) +
    scale_y_continuous(expand = c(0, 0), breaks = 2) +
    ggtitle(title, paste0(subtitle, "\n", 
                          paste(unique(dat[!is.na(dat$n_degs), ]$n_degs), collapse = ", "))) +
    labs(y = "-log10(adjusted p-value)", x = "logFC")
}

# pseudobulks ####

# load pseudobulks
pse <- list()
system(
  "ls output/snRNAseq_workflow/by_patient_wo_organoids/*/pseudobulk/pse.rds", 
  intern = T) %>%
  purrr::walk(function(file) {
    patient <- stringr::str_split(file, "/")[[1]][4]
    print(patient)
    pse[[patient]] <<- readRDS(file)
  })

# combine pseudobulks
pse <- mergeSEs(pse, do.scale = F)

# transform some variables
pse$broad_response <- ifelse(
  pse$response %in% c("Responder", "Partial Responder"),
  "responder", "nonresponder"
)
pse$malignancy_quartile <-
  cut(pse$malignancy, quantile(pse$malignancy),
      include.lowest = T, labels = F)
pse$cycling_quartile <-
  cut(pse$cycling_score, quantile(pse$cycling_score),
      include.lowest = T, labels = F)
pse$lineage <- ifelse(
  pse$label.infercnv == "Malignant_cells", 
  "malignant",
  ifelse(
    pse$label.infercnv %in% 
      c("NK_cells", "Myeloid_cells", "Erythrocytes", 
        "Macrophage", "T_cells", "B_cells"), 
    "immune",
    "stromal")
)

# save
dir.create("output/glmGamPoi/")
saveRDS(pse, "output/glmGamPoi/pse.rds")

# glmGamPoi test_de
p <- list()
degs <- list()

# fit model
fit <- glm_gp(
  pse, 
  design = ~ 
    label.infercnv + 
    timepoint + 
    broad_response + 
    sample_site +
    malignancy_quartile + 
    cycling_quartile + 
    age + 
    gender 
  )

# DE of nonresponsive malignant cells (vs responsive)
degs[["nonresponse"]] <- test_de(
  fit, contrast = 
    cond(label.infercnv = "Malignant_cells", broad_response = "nonresponder") -
    cond(label.infercnv = "Malignant_cells", broad_response = "responder"), 
  sort_by = adj_pval 
)

# DE of treated malignant cells (vs untreated)
degs[["treatment"]] <- test_de(
  fit, contrast = 
    cond(label.infercnv = "Malignant_cells", timepoint = "post") -
    cond(label.infercnv = "Malignant_cells", timepoint = "pre"), 
  sort_by = pval) 

# DE of metastasised malignant cells (vs local)
degs[["metastasis"]] <- test_de(
  fit, contrast = 
    cond(label.infercnv = "Malignant_cells", sample_site = "lymph_node") -
    cond(label.infercnv = "Malignant_cells", sample_site = "tumour"), 
  sort_by = adj_pval)

# DE of high-malignancy malignant cells (vs low-malignancy)
degs[["malignancy"]] <- test_de(
  fit, contrast = 
    cond(label.infercnv = "Malignant_cells", malignancy_quartile = 4) -
    cond(label.infercnv = "Malignant_cells", malignancy_quartile = 1),
  sort_by = adj_pval)

# DE of malignant cells (vs epithelial)
degs[["transformation"]] <- test_de(
  fit, contrast = 
    cond(label.infercnv = "Malignant_cells") -
    cond(label.infercnv = "Epithelial_cells"),
  sort_by = adj_pval
)

# DE of male malignant cells (vs female)
degs[["sex"]] <- test_de(
  fit, contrast = 
    cond(label.infercnv = "Malignant_cells", gender = "M") -
    cond(label.infercnv = "Malignant_cells", gender = "F"),
  sort_by = adj_pval
)

# plots
p[["nonresponse"]] <-
  plot_volcano(degs$nonresponse, "DEGs of treatment resistance",
               "treatment-resistant malignant cells (vs responsive)")
p[["treatment"]] <-
  plot_volcano(degs$treatment, "DEGs of treatment", 
               "treated malignant cells (vs untreated)")
p[["metastasis"]] <-
  plot_volcano(degs$metastasis, "DEGs of metastasis", 
               "lymph node malignant cells (vs primary)")
p[["malignancy"]] <-
  plot_volcano(degs$malignancy, "DEGs of malignancy",
               "high-malignancy malignant cells (vs low)")
p[["sex"]] <-
  plot_volcano(degs$sex, "DEGs of sex",
               "male malignant cells (vs female)")

# save plots
p %>%
  names() %>%
  purrr::walk(function(test) {
    pdf(paste0("output/glmGamPoi/", test, "_volcano.pdf"), width = 5)
    print(p[[test]])
    dev.off()
  })

# save degs
saveRDS(degs, "output/glmGamPoi/degs.rds")
degs %>% 
  names() %>%
  purrr::walk(function(test) {
    readr::write_tsv(degs[[test]], paste0("output/glmGamPoi/", test, "_degs.tsv"))
    degs[[test]] %>% 
      tibble() %>%
      transmute(test = test,
                gene = name, adj_pval, logFC = lfc, 
                direction = ifelse(logFC < 0, "downregulated", "upregulated")) %>%
      filter(abs(logFC) > 1.5, adj_pval < 0.01) %>%
      group_by(direction) %>%
      slice_max(logFC, n = 10) %>%
      readr::write_tsv(paste0("output/glmGamPoi/", test, "_top_degs.tsv"))
  })
degs %>%
  bind_rows(.id = "test") %>%
  transmute(test = test,
            gene = name, adj_pval, logFC = lfc, 
            direction = ifelse(logFC < 0, "downregulated", "upregulated")) %>%
  filter(abs(logFC) > 1.5, adj_pval < 0.01) %>%
  group_by(test, direction) %>%
  slice_max(logFC, n = 10) %>%
  readr::write_tsv("output/glmGamPoi/top_degs_by_logFC.tsv")
degs %>%
  bind_rows(.id = "test") %>%
  transmute(test = test,
            gene = name, adj_pval, logFC = lfc, 
            direction = ifelse(logFC < 0, "downregulated", "upregulated")) %>%
  filter(abs(logFC) > 1.5, adj_pval < 0.01) %>%
  group_by(test, direction) %>%
  slice_min(adj_pval, n = 10) %>%
  readr::write_tsv("output/glmGamPoi/top_degs_by_adj_pval.tsv")
  
  
