base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "/rds/general/user/art4017/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

dir.create("output/celltype_composition/")

library(magrittr)
library(ggplot2)

# celltype colours
colours <- c(
  "Malignant_cells" = "black",
  "Endothelial_cells" = "#AD7700",
  "Epithelial_cells" = "#ffdcc4",
  "Fibroblasts" = "#CC79A7",
  "Smooth_muscle_cells" = "#fae64d",
  "Tissue_stem_cells" = "#D55E00",
  "Stromal_cells" = "#ffc02e",
  # immune
  "NK_cells" = "#0072B2",
  "T_cells" = "#4ea8ad",
  "B_cells" = "#009E73",
  "Macrophage" = "#00F6B3",
  "Myeloid_cells" = "#ce95fc",
  "Erythrocytes" = "#d45353",
  # other
  "Others" = "#bfbfbf"
)

# infercnv annots
ref_annots <- c("T_cells", "Macrophage", "B_cells", "NK_cells", "Myeloid_cells", "Stromal_cells")
query_annots <- "Epithelial_cells"
ref_and_query <- 
  list("ref" = ref_annots, "query" = query_annots) %>% 
  tibble::enframe(name = "group", value = "label.reduced") %>%
  tidyr::unnest(cols = "label.reduced")

# before infercnv ####

# # read in singler annots
# singler_annots_files <-
#   system(
#     "ls output/snRNAseq_workflow/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds",
#     intern = T
#   )
# singler_annots <- list()
# singler_annots_files %>%
#   purrr::walk(function(file) {
#     patient <- strsplit(file, "/")[[1]][4]
#     print(patient)
#     singler_annots[[patient]] <<- readRDS(file)@meta.data
#   })
# singler_annots <-
#   singler_annots %>%
#   dplyr::bind_rows(.id = "patient") %>%
#   tibble::as_tibble() %>%
#   dplyr::mutate(
#     label.reduced = factor(label.reduced, levels = names(colours)),
#     timepoint = factor(timepoint, levels = c("pre", "post"))
#   )
# 
# # save data
# readr::write_tsv(singler_annots, "output/celltype_composition/singler_annots.tsv")
# 
# # barplot
# pdf("output/celltype_composition/singler_composition.pdf", height = 7, width = 15, onefile = T)
# singler_annots %>%
#   ggplot(aes(x = sample, fill = label.reduced)) +
#   geom_bar(position = "fill") +
#   scale_fill_manual(values = colours) +
#   facet_grid(~ patient, scales = "free", space = "free_x") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
#         axis.text.y = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_blank())
# singler_annots %>%
#   ggplot(aes(x = sample, fill = label.reduced)) +
#   geom_bar() +
#   scale_fill_manual(values = colours) +
#   facet_grid(~ patient, scales = "free", space = "free_x") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
#         axis.text.y = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_blank())
# dev.off()

# after infercnv ####

# read in singler annots
singler_annots_files <-
  system(
    "ls output/snRNAseq_workflow/by_patient_wo_organoids/*/integrating/infercnv/seu.rds",
    intern = T
  )
singler_annots <- list()
singler_annots_files %>%
  purrr::walk(function(file) {
    patient <- strsplit(file, "/")[[1]][4]
    print(patient)
    singler_annots[[patient]] <<- readRDS(file)@meta.data
  })
singler_annots <-
  singler_annots %>%
  dplyr::bind_rows(.id = "patient") %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    label.infercnv = factor(label.infercnv, levels = names(colours)),
    timepoint = factor(timepoint, levels = c("pre", "post")),
    lineage = dplyr::case_when(
      label.infercnv == "Malignant_cells" ~ "malignant",
      label.infercnv %in% c("NK_cells", "Myeloid_cells", "Erythrocytes", 
                            "Macrophage", "T_cells", "B_cells") ~ "immune",
      .default = "stromal"
    )
  )

# save data
readr::write_tsv(singler_annots, "output/celltype_composition/infercnv_singler_annots.tsv")

# barplot
pdf("output/celltype_composition/infercnv_singler_composition.pdf", height = 7, width = 15, onefile = T)
singler_annots %>%
  ggplot(aes(x = sample, fill = label.infercnv)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = colours) +
  facet_grid(~ patient, scales = "free", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
singler_annots %>%
  ggplot(aes(x = sample, fill = label.infercnv)) +
  geom_bar() +
  scale_fill_manual(values = colours) +
  facet_grid(~ patient, scales = "free", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
singler_annots %>%
  ggplot(aes(x = gender, fill = label.infercnv)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = colours) +
  facet_grid(~ sample_site, scales = "free", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
singler_annots %>%
  ggplot(aes(x = gender, fill = label.infercnv)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = colours) +
  facet_grid(lineage ~ sample_site, scales = "free", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
dev.off()

