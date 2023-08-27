patient <- "A"

# get integrated seu
seu <- 
  readRDS(
    paste0("output/snRNAseq_workflow/by_patient_wo_organoids/", 
           patient, "/integrating/infercnv/malignant/clustering/seu_annotated.rds")
  )

# convert to cds
Seurat::DimPlot(seu, group.by = "SCT_snn_res.0.1")

# find markers of pre- vs post-treatment
tests <- list(timepoint = "post", sample_site = "lymph_node", label.infercnv = "Malignant_cells")
markers <-
  tests %>%
  names() %>%
  purrr::map(
    function(test) {
      print(test)
      seu %>%
        Seurat::FindMarkers(
          ident.1 = colnames(seu[, seu@meta.data[, test] == tests[[test]]]), 
          test.use = "MAST"
        ) %>%
        tibble::as_tibble(rownames = "gene")
    }
  ) %>%
  setNames(unlist(tests)) %>%
  dplyr::bind_rows(.id = "test")
markers %>%
  dplyr::arrange(-abs(avg_log2FC))

malig_vs_epith <-
  seu[, seu$label.infercnv %in% c("Epithelial_cells", "Malignant_cells")] %>%
  Seurat::FindMarkers(
    ident.1 = colnames(seu[, seu$label.infercnv == "Malignant_cells"]),
    test.use = "MAST"
  )
malig_vs_epith %>%
  tibble::as_tibble(rownames = "gene") %>%
  dplyr::arrange(-abs(avg_log2FC))
