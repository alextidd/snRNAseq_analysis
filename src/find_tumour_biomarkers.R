base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "/rds/general/user/art4017/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(ggVennDiagram)

plot_markers <- function(markers, top_n = 10) {
  p_dat <-
    markers %>%
    group_by(direction, n) %>%
    mutate(rank = rank(-abs(avg_log2FC))) %>%
    filter(rank <= top_n) %>%
    mutate(`-log10(p-value)` = -log10(p_val_adj_no_0))
  p_dat %>%
    ggplot(aes(x = -rank, y = avg_log2FC, 
               fill = direction, alpha = `-log10(p-value)`,
               label = gene)) +
    geom_bar(stat = "identity") +
    geom_text(hjust = ifelse(p_dat$group == 2, 1.1, -0.1),
              alpha = 1,
              size = 4) +
    scale_y_continuous(labels = abs, limits = max(abs(p_dat$avg_log2FC)) * c(-1.3, 1.3)) +
    scale_fill_manual(values = as.vector(c("#3179de", "#d23f67"))) +
    scale_alpha(range = c(0.5, 1)) +
    labs(x = "", y = "", fill="") +
    theme_minimal() +   
    coord_flip() +
    theme( 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size=10), 
      axis.text.y = element_blank(),
      strip.text.x=element_text(size=15),
      legend.position="bottom",
      legend.text=element_text(size=10)
    ) 
}
fix_pathway_names <- function(df) {
  df %>% 
    mutate(Pathway = Pathway %>% tolower() %>% gsub("\\_", " ", .) %>%
             gsub("hallmark", "H", .) %>%
             gsub("reactome", "R", .) %>%
             gsub("KEGG", "K", .))
}
plot_my_heatmap <- function(p_dat_in, col, stat = "qval") {
  sig <- 
    p_dat_in %>%
    select(Pathway, value = adjPval, name = matches(col)) %>%
    pivot_wider() %>%
    column_to_rownames("Pathway") %>%
    as.matrix() 
  p_dat <-
    p_dat_in %>%
    select(Pathway, value = all_of(stat), name = all_of(col)) %>%
    pivot_wider() %>%
    column_to_rownames("Pathway") %>%
    as.matrix() 
  Heatmap(
    p_dat,
    name = stat,
    border = T,
    show_row_dend = F,
    show_column_dend = T,
    row_names_gp = grid::gpar(fontsize = 8), 
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb = textGrob("*")
      gb_w = convertWidth(grobWidth(gb), "mm")
      gb_h = convertHeight(grobHeight(gb), "mm")
      if (sig[i, j] < 0.0005) {
        grid.text("***", x, y - gb_h*0.5 + gb_w*0.4)
      } else if(sig[i, j] < 0.005) {
        grid.text("**", x, y - gb_h*0.5 + gb_w*0.4)
      } else if (sig[i, j] < 0.05) {
        grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
      }
    }     
  )
}
plot_my_enrichment <- function(p_dat_in, top_n = 5) {
  p_dat <-
    p_dat_in %>%
    group_by(patient) %>%
    mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                             FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                             FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                             FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'),
           rank = rank(-qval, ties.method = "min"),
           lab = case_when(
             rank < 4 & qval > 0 & abs(FC) > 5 & adjPval < 0.01 ~ Pathway,
             TRUE ~ ""))
  
  ggplot(p_dat, aes(FC, qval)) +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
    geom_point(cex = 2.6, shape = 21, fill = p_dat$color, stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = lab),
                             size = 4, min.segment.length = 0, nudge_x = 10) +
    ylab("Qval") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          aspect.ratio = 1) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
}

# drivers
drivers <- c("TP53", "CDKN2A", "SMAD4", "ARID1A", "ERBB2", "KRAS", "PIK3CA", "SMARCA4", "CTNNB1", "ARID2", "PBRM1", "FBXW7")
oac_markers <- c("EPCAM", "AGR2", "CEACAM6", "KRT8", "KRT18", "KRT19")

# read in patient metadata
sample_metadata <-
  readr::read_tsv("../cellranger/output/sample_metadata.tsv") 
patient_metadata <-
  sample_metadata %>%
  transmute(
    patient = patient_id,
    broad_response = case_when(response %in% c("Responder", "Partial Responder") ~ "Responder",
                               response %in% c("Progressor", "Non-responder") ~ "Non-responder"),
    response,
    gender,
    age,
    treatment) %>%
  distinct()

# read in degs
degs <- list()
system("ls output/snRNAseq_workflow/by_patient_wo_organoids/*/integrating/differential_expression/degs.rds",
         intern = T) %>%
  purrr::walk(function(path) {
    patient <- strsplit(path, "/")[[1]][4]
    degs[[patient]] <<- readRDS(path)
  })
degs <-
  degs %>% 
  purrr::map(~ bind_rows(.x, .id = "test")) %>%
  bind_rows(.id = "patient") %>%
  left_join(patient_metadata)

# read in pathways
scpa <- list()
system("ls output/snRNAseq_workflow/by_patient_wo_organoids/*/integrating/differential_expression/differential_expression_cache/scpa_*.rds",
       intern = T) %>%
  purrr::walk(function(path) {
    patient <- strsplit(path, "/")[[1]][4]
    scpa[[patient]] <<- readRDS(path)$malignant_vs_epithelial
  })
scpa <-
  scpa %>%
  bind_rows(.id = "patient") %>%
  fix_pathway_names()

# read in malignant cluster markers
malig_clust_markers <-
  readRDS("output/snRNAseq_workflow/by_patient_wo_organoids/malignant/clustering/clustering_cache/cluster_markers_a6b5f1f8c7b8f832913507e49f1b1038.rds")

# top 50 pathways
top50 <-
  scpa %>%
  group_by(Pathway) %>%
  summarise(qval = max(qval)) %>%
  slice_max(order_by = qval, n = 50) %>%
  pull(Pathway)
most_common <-
  scpa %>%
  filter(adjPval < 0.01) %>%
  group_by(Pathway) %>%
  summarise(n = n_distinct(patient)) %>%
  arrange(-n) %>%
  slice_max(order_by = n, n = 50) %>%
  pull(Pathway)
hallmark <- 
  scpa %>% 
  filter(grepl("H ", Pathway))
plot_my_heatmap(scpa %>% filter(Pathway %in% most_common), "patient")
plot_my_heatmap(hallmark, "patient")

# plot heatmap
dir.create("output/pathway_analysis/")
pdf("output/pathway_analysis/scpa_malignant_vs_epithelial_heatmap.pdf", width = 10)
gridExtra::marrangeGrob(
  list(grid.grabExpr(draw(plot_my_heatmap(scpa, "patient"))),
       grid.grabExpr(draw(plot_my_heatmap(scpa, "patient", stat = "FC")))),
  nrow = 1, ncol = 2, top = NULL
)
dev.off()

# plot enrichment
pdf("output/pathway_analysis/scpa_malignant_vs_epithelial.pdf", 
    width = 12, height = 10)
scpa %>%
  plot_my_enrichment() + 
  facet_wrap(~patient) 
dev.off()

# plot shared biomarkers across patients
pdf("output/pathway_analysis/shared_degs.pdf",
    width = 7,
    height = 12)
degs %>%
  filter(test == "malignant_vs_epithelial") %>%
  group_by(gene, direction) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  summarise(across(where(is.numeric), mean)) %>%
  arrange(-n,-abs(avg_log2FC)) %>%
  mutate(n = paste(n, "patients")) %>%
  plot_markers() +
  ggforce::facet_col( ~ n,
                      scales = "free_y",
                      space = "free") +
  theme(strip.text.x = element_text(size = 12))
dev.off()

# plot venn diagram of marker sharing
# degs %>% 
#   group_by(gene, direction) %>%
#   filter(test == "malignant_vs_epithelial",
#          direction == "Malignant_cells") %>%
#   {split(.$gene, .$broad_response)} %>%
#   ggVennDiagram()

# plot malignant cluster markers
pdf("output/pathway_analysis/malignant_cluster_markers.pdf",
    width = 15, height = 15)
malig_clust_markers %>%
  tibble() %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  mutate(cluster_lab = paste0("cluster ", cluster, " (", n(), " markers)")) %>%
  {mutate(., cluster_lab = factor(
    cluster_lab, 
    levels = unique(arrange(., cluster)$cluster_lab)))} %>%
  mutate(lab = case_when(rank(-avg_log2FC) <= 5 ~ gene,
                         TRUE ~ "")
         ) %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj),
             label = lab)) +
  geom_point(alpha = 0.5) +
  ggrepel::geom_text_repel(max.overlaps = 20) +
  facet_wrap(~cluster_lab) +
  theme_classic() +
  theme(strip.text = element_text(size = 20))
dev.off()

# plot shared markers
degs %>%
  filter(test == "malignant_vs_epithelial") %>%
  group_by(gene, direction) %>% 
  mutate(n = n()) %>%
  filter(n > 4) %>%
  arrange(direction, desc(n)) %>%
  ggplot(aes(x = patient, y = reorder(gene, direction), fill = direction)) +
  geom_tile() +
  theme_classic() +
  facet_wrap(~ direction, scales = "free_y")

p_dat <-
  degs %>%
  filter(abs(avg_log2FC) > 0.75) %>%
  group_by(gene) %>%
  arrange(n()) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  distinct(gene, value = avg_log2FC, name = patient) %>%
  pivot_wider() %>%
  column_to_rownames("gene") %>%
  as.matrix()
p_ann <-
  tibble(patient = colnames(p_dat)) %>%
  left_join(patient_metadata) 
p_dat  %>%
  Heatmap(
    row_names_gp = grid::gpar(fontsize = 6),
    top_annotation = HeatmapAnnotation(
      response = p_ann$broad_response,
      treatment = p_ann$treatment
    )
  )

degs %>%
  filter(abs(avg_log2FC) > 0.75) %>%
  group_by(gene) %>%
  filter(all(c("Malignant_cells", "Epithelial_cells") %in% direction),
         n() > 1) %>%
  group_by(gene, direction) %>%
  summarise(patient = paste(patient, collapse = ","),
            p_val = mean(p_val),
            avg_log2FC = mean(avg_log2FC))


p_dat <- degs %>%
  filter(abs(avg_log2FC) > 0.75) %>%
  group_by(gene) %>%
  filter(all(c("Malignant_cells", "Epithelial_cells") %in% direction)) %>%
  select(gene, value = avg_log2FC, name = patient) %>%
  pivot_wider(values_fill = 0) %>%
  column_to_rownames("gene") %>%
  as.matrix() 
p_ann <-
  tibble(patient = colnames(p_dat)) %>%
  left_join(patient_metadata) 
p_dat %>%
  Heatmap(top_annotation = HeatmapAnnotation(broad_response = p_ann$broad_response))

# genes with divergent differential expression direction between different patients
degs %>%
  #filter(abs(avg_log2FC) > 0.75) %>%
  group_by(gene) %>%
  filter(all(c("Malignant_cells", "Epithelial_cells") %in% direction)) %>%
  group_by(gene, direction) %>%
  summarise(avg_log2FC = mean(avg_log2FC), 
            p_val_adj = mean(p_val_adj),
            patient = paste(patient, collapse = ","),
            broad_response = paste(unique(broad_response), collapse = ",")) %>%
  group_by(gene) %>% 
  mutate(differential = abs(max(avg_log2FC) - min(avg_log2FC))) %>%
  arrange(-differential) %>%
  readr::write_tsv("output/pathway_analysis/divergent_degs.tsv")

