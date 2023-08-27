patient <- "A"

# get integrated seu
seu <- readRDS(
  paste0("output/snRNAseq_workflow/by_patient_wo_organoids/", 
         patient, "/integrating/seurat_clustering/seu_annotated.rds")
  )

# get aneuploid cells
input_dir <- paste0(
  "output/snRNAseq_workflow/by_patient_wo_organoids/", patient, 
  "/integrating/infercnv/infercnv_cache/"
)
samples <- list.files(input_dir)
malignant_cells <-
  samples %>%
  purrr::map(function(sample_i) {
    readr::read_tsv(
      paste0(input_dir, sample_i, 
             "/HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")) %>%
      dplyr::mutate(sample = sample_i)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(sample, cell_group_name) %>%
  dplyr::summarise(aneuploid = any(c(1, 3) %in% state)) %>%
  dplyr::filter(
    aneuploid & cell_group_name %in% colnames(seu[, seu$label.reduced == "Epithelial_cells"])
  ) %>%
  dplyr::pull(cell_group_name)

# add malignant label
seu$label.infercnv <- seu$label.reduced
seu@meta.data[malignant_cells, "label.infercnv"] <- "Malignant_cells"

# propeller
prop_run <-
  seu[, seu$ == "tumour"] %>% {
    speckle::propeller(
      clusters = .$label.infercnv,
      sample = .$sample,
      group = .$sample_site,
      robust = F, trend = F, transform = "asin"
    )
  }
  

  
