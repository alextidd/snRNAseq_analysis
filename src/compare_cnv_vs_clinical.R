base_dir <- ifelse(Sys.info()["login"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "~/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

# load sample metadata
sample_metadata <-
  readr::read_tsv("../cellranger/output/sample_metadata.tsv")

# load predictions
cnv_files <- 
  list.files(
    "output/snRNAseq_workflow/by_patient_wo_organoids/",
    pattern = "infercnv/seu_infercnv_malignant.rds",
    recursive = T
  ) %>%
  setNames(gsub("/.*", "", .))
cnv_preds <- 
  cnv_files %>% 
  names() %>%
  purrr::map(function(sample_i) {
    print(sample_i)
    readRDS(
      paste0("output/snRNAseq_workflow/by_patient_wo_organoids/", 
             cnv_files[[sample_i]]))@meta.data
  }) %>%
  dplyr::bind_rows() 

seu_infercnv@meta.data %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(label.infercnv == "Malignant_cells") %>% 
  dplyr::select(dplyr::starts_with("has_")) %>% 
  tidyr::pivot_longer(everything()) %>% 
  dplyr::mutate(value = as.numeric(value)) %>% 
  dplyr::filter(grepl("dupli|loss", name)) %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(prop = 100 * sum(value) / dplyr::n()) %>%
  dplyr::arrange(-prop) %>% 
  dplyr::mutate(
    chr = gsub(".*chr", "", name) %>% as.numeric(),
    direction = gsub("has_", "", name) %>% gsub("_.*", "", .),
    y = ifelse(direction == "dupli", 1, -1) * prop
  )  %>% 
  ggplot2::ggplot(ggplot2::aes(x = chr,  y = y, fill = direction)) + 
  ggplot2::geom_col() +
  ggplot2::labs(y = "percentage gain or loss")

