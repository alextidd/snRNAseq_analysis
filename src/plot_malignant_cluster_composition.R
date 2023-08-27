base_dir <- ifelse(Sys.info()["user"] == "alexandratidd1",
                   "/Volumes/art4017/",
                   "/rds/general/user/art4017/")
setwd(paste0(base_dir,"home/snRNAseq_analysis/"))

dir.create("output/malignant_clusters/")

library(magrittr)
library(dplyr)
library(ggplot2)

cluster_metadata <-
  readRDS("output/snRNAseq_workflow/by_patient_wo_organoids/malignant/clustering/clustering_cache/cluster_metadata_9869b779e7071db7d5e9a232bcca5aec.rds")

pdf("output/malignant_clusters/patient_composition.pdf", 
    height = 5)
cluster_metadata %>%
  ggplot(aes(x = cluster, fill = patient_id)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  theme_classic()
dev.off()

cluster_markers <-
  readRDS("output/snRNAseq_workflow/by_patient_wo_organoids/malignant/integrating_mnn/integrating_mnn_cache/markers_aca116f7396ff51fccab829d8568d13d.rds")

