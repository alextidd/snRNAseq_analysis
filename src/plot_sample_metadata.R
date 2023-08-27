base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "/rds/general/user/art4017/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

library(magrittr)

# get sample metadata
sample_metadata <-
  readr::read_tsv(paste0(base_dir, "/home/cellranger/output/sample_metadata.tsv")) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(n_samples = dplyr::n(), 
                treatment = ifelse(treatment == "NONE", "unknown", treatment)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    readr::read_tsv("output/celltype_composition/singler_annots.tsv") %>%
      dplyr::count(sample, name = "n_cells")
  )

# response colours / levels
response_colours <-
  c("Responder" = "lightgreen",
    "Partial Responder" = "#ffe97d",
    "Non-responder" = "#ffb24d",
    "Progressor" = "#ff9c9c")

# heatmap
p_dat <-
  sample_metadata %>%
  dplyr::filter(sample_type != "organoid") %>%
  dplyr::transmute(
    # patient-level variables
    patient_id,
    response = factor(response, levels = names(response_colours)),
    gender, age, n_samples,
    treatment,
    # sample-level variables
    sample, timepoint, Tx, Mx, Nx, sample_site, mandard_score, batch_dir, n_cells
  ) %>%
  dplyr::arrange(response, patient_id, desc(timepoint))
p_dat$sample <- factor(p_dat$sample, levels = unique(p_dat$sample))
p_dat$patient_id <- factor(p_dat$patient_id, levels = unique(p_dat$patient_id))

plot_variable <- function(p_dat, variable, lvl, variable_colours,
                          show_col_names = F, return_legend = F) {
  p <-
    p_dat %>%
    dplyr::select(sample, dplyr::all_of(variable)) %>%
    dplyr::mutate(facet = p_dat$patient_id) %>%
    tidyr::pivot_longer(-c(sample, facet)) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = sample, y = name, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::theme_void() +
    ggplot2::facet_grid(~facet, scales = "free_x", space = "free_x", ) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1, size = 15), 
                   strip.text.x = ggplot2::element_blank(),
                   panel.spacing.x = ggplot2::unit(0.5, "points"))
  if (show_col_names == T) {
    p <- p +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 15)) 
  }
  if (all(variable %in% names(variable_colours))) {
    p <- p +
      ggplot2::scale_fill_manual(values = variable_colours[[variable]],
                                 na.value = "white")
  } else if (is.numeric(unlist(p_dat[, variable , drop=T]))) {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#fae3f5",
      high = "#611d52",
      na.value = 'white')
  } else {
    p <- p + ditto_colours
  }
  if (return_legend == T) {
    lemon::g_legend(
      p +
        ggplot2::labs(fill = ifelse(lvl == "mut", "mut ccf", variable)) +
        ggplot2::theme(legend.position = "bottom")
    )
  } else {
    p + ggplot2::theme(legend.position = "none") 
  }
}

# variable levels
variable_lvls <- 
  list("patient_id", "gender", "age", "batch_dir", 
       "timepoint", "sample_site", "Tx", "Nx", "n_cells",
       "treatment", "response")

# get colours
ditto_colours <- list(ggplot2::scale_fill_manual(values = dittoSeq::dittoColors(), na.value = "white"),
                      ggplot2::scale_colour_manual(values = dittoSeq::dittoColors(), na.value = "white"))
avail_colours <- dittoSeq::dittoColors()
variable_colours <- list()
p_leg <- list()
for(variable in variable_lvls) {
  if (!is.numeric(unlist(p_dat[,variable]))) {
    cat(variable, "categorical\n")
    n_colours <- length(unique(p_dat[,variable,drop=T]))
    variable_colours[[variable]] <-
      avail_colours[1:n_colours]
    avail_colours <- avail_colours[-c(1:n_colours)]
  }
}
variable_colours[["gender"]] <- c("F" = "#ffb3fc", "M" = "#7eaaed")
variable_colours[["response"]] <- response_colours

# patient / sample level variables
p <- list()
p[["patient"]] <- 
  plot_variable(p_dat, "patient_id", "sample", variable_colours, show_col_names = T)
variable_lvls[-1] %>% unlist() %>% 
  purrr::walk(function(variable) {
    cat(variable, "\n")
    p_leg[[variable]] <<-
      plot_variable(p_dat, variable, "sample", variable_colours, return_legend = T)
    p[[variable]] <<-
      plot_variable(p_dat, variable, "sample", variable_colours)
  })

# save plots
pdf(paste0("output/sample_heatmap.pdf"), width = 12, height = 3.5, onefile = T)
patchwork::wrap_plots(p, ncol = 1)
dev.off()

pdf(paste0("output/sample_heatmap_legends.pdf"), width = 6, height = 10)
p_leg %>% patchwork::wrap_plots(ncol = 1)
dev.off()