markers <- 
  readr::read_tsv("../sc_studies/data/markers/32532891_2021_Zhang/markers.tsv")
avail_markers <- 
  markers %>%
  dplyr::filter(gene %in% rownames(seu)) %>%
  {split(.$gene, .$population)}
reduced_markers <- list(
  "T_cells" = avail_markers$`T cells`,
  "B_cells" = avail_markers$`B cells`,
  "Macrophage" = avail_markers$macrophages,
  "Epithelial_cells" = avail_markers$epithelium,
  "Fibroblasts" = avail_markers$fibroblasts
)

dittoSeq::dittoDimPlot(seu, "label.infercnv", cells.use = colnames(seu[, seu$label.infercnv == "B_cells"])) +
  Nebulosa::plot_density(seu, avail_markers[[ct]], joint = T, combine = F) %>%
  tail(1) 
