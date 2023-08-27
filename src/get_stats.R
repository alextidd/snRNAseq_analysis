
library(ggplot2)
library(tidyverse)
library(magrittr)
library(tidytext)

base_dir <- ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/art4017/home/","~/")
setwd(paste0(base_dir, "snRNAseq_analysis/"))

# read
sample_metadata <- read_tsv(paste0(base_dir,"cellranger/output/sample_metadata.tsv"))
metrics <-
  sample_metadata$sample %>%
  unique() %>%
  map(function(samp) {
    metrics_file <- paste0(base_dir, "cellranger/output/", samp, "/outs/metrics_summary.csv")
    if (file.exists(metrics_file)) {
      metrics_file %>%
        read_csv(show_col_types = F) %>%
        mutate(sample = samp) %>%
        select(sample, everything())
    } else { tibble() }
  }) %>%
  bind_rows() %>%
  left_join(sample_metadata) %>%
  mutate(sample_code = factor(sample_code, levels = c(paste0("T", 0:3), paste0("N", 0:3))),
         timepoint = factor(timepoint, levels = c("pre", "post")))

statistics <- 
  c("Estimated Number of Cells", "Mean Reads per Cell", "Median Genes per Cell")

pdf("output/n_nuclei.pdf", width = 16,height = 8, onefile = T)
metrics %>% 
  pivot_longer(cols = statistics,
               names_to = "statistic") %>%
  ggplot(aes(x = sample_code, y = value, fill = timepoint)) +
  geom_col() + 
  theme_classic() +
  facet_grid(statistic ~ patient_id, scales = "free", space = "free_x")  +
  scale_fill_brewer(palette = "Dark2")
metrics %>%
  group_by(patient_id) %>%
  mutate(total_n_nuclei = sum(`Estimated Number of Cells`),
         label = prettyNum(total_n_nuclei, big.mark=",", scientific=FALSE)) %>%
  ggplot(aes(x = patient_id, y = `Estimated Number of Cells`, 
             fill = sample_code, label = label)) +
  geom_col(position = "stack") + 
  geom_label(aes(y = total_n_nuclei), fill = "white") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")
metrics %>%
  pivot_longer(cols = c("Median Genes per Cell", "Mean Reads per Cell"),
               names_to = "statistic") %>%
  ggplot(aes(x = sample, y = value)) +
  geom_col() + 
  theme_classic() +
  facet_grid(statistic ~ sample_type, scales = "free", space = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 14)) 
dev.off()

metrics %>%
  select(sample, batch_dir, all_of(statistics)) %>%
  write_tsv("output/metrics.tsv")



