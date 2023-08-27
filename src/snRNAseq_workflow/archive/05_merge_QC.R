setwd("~/snRNAseq_analysis/")

data_dir <- "output/snRNAseq_workflow/by_sample/"

# read in qc
samples <-
  list.files(path = data_dir, include.dirs = T)
fail_criteria_stats <-
  samples %>%
  purrr::map(function(sample) {
    filts <-
      paste0(data_dir, "/", sample, "/filtering/filters.rds") %>%
      readRDS()
    n_retained <- nrow(dplyr::filter(filts$nucleus_filtering, pass))
    n_total <- nrow(filts$nucleus_filtering)
    dplyr::left_join(
      filts$filters %>% dplyr::bind_rows(.id = "fail_criteria") %>%
        tidyr::pivot_longer(c("min", "max")) %>%
        dplyr::transmute(fail_criteria = paste(name, fail_criteria, sep = "_"), value),
      filts$nucleus_filtering %>%
        dplyr::filter(!pass) %>%
        tidyr::separate_longer_delim(fail_criteria, ",") %>%
        dplyr::count(fail_criteria) %>%
        dplyr::mutate(`%` = round((n / n_total) * 100, 2)),
      by = "fail_criteria"
    )
  }) %>% 
  setNames(samples) %>%
  dplyr::bind_rows(.id = "sample")

# nucleus filtering
nucleus_filtering <-
  samples %>%
  purrr::map(function(sample) {
    filts <-
      paste0(data_dir, "/", sample, "/filtering/filters.rds") %>%
      readRDS()
    filts$nucleus_filtering
  }) %>% 
  setNames(samples) %>%
  dplyr::bind_rows()

# plot
p_dat <-
  nucleus_filtering %>%
  dplyr::group_by(sample) %>%
  dplyr::count(fail_criteria) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::left_join(
    nucleus_filtering %>%
      dplyr::filter(is.na(fail_criteria)) %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(pass = dplyr::n())
  ) %>%
  dplyr::mutate(
    patient = sample %>% gsub("\\_.*", "", .),
    label = paste0(pass, " / ", total, " passed (", round((pass/total)*100, 1), "%)")
  )
pdf("output/snRNAseq_workflow/by_sample/merge_qc.pdf", width = 20, height = 20)
p_dat %>%
  head(10000) %>%
  ggplot2::ggplot(ggplot2::aes(x = sample, y = n, fill = fail_criteria, label = label)) +
  ggplot2::geom_col() +
  ggplot2::geom_label(ggplot2::aes(y = total + 500), fill = "white", hjust = 0) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 7)) +
  ggplot2::scale_fill_manual(values = dittoSeq::dittoColors(), na.value = "grey") +
  ggplot2::coord_flip() +
  ggplot2::ylim(0, max(p_dat$total) + 13000) +
  ggplot2::facet_grid(rows = "patient", scales = "free_y", space = "free_y")
dev.off()

