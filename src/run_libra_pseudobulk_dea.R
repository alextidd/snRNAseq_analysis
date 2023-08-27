library(magrittr)

seu <- readRDS("output/snRNAseq_workflow/by_patient_wo_organoids/A/integrating_mnn/seu.rds")
meta <- readRDS("output/snRNAseq_workflow/by_patient_wo_organoids/A/integrating/infercnv/seu.rds")@meta.data
seu$label <- paste(seu$patient_id, seu$sample_site, sep = "_")
seu$replicate <- seu$sample_code
seu$cell_type <- meta$label.infercnv

Seurat::DefaultAssay(seu) <- "RNA"
de <-
  Libra::run_de(seu)
    
de %>%
  dplyr::filter(p_val_adj < 0.05, avg_logFC > 0) %>%
  dplyr::arrange(-avg_logFC)

de_malig_vs_epith <-
  Libra::run_de(seu[, seu$cell_type %in% c("Epithelial_cells", "Malignant_cells")])
de_malig_vs_epith %>%
  dplyr::filter(p_val_adj < 0.05, avg_logFC > 0) %>%
  dplyr::arrange(-avg_logFC)

seu$label <- seu$timepoint
de_pre_vs_post <-
  Libra::run_de(seu[, seu$cell_type == "Malignant_cells"])

# pseudobulk with AverageExpression
Seurat::Idents(seu) <- seu$label.infercnv
seu <- Seurat::ScaleData(seu)
mat <- Seurat::AverageExpression(
  seu, group.by = c("sample_code", "timepoint"))


