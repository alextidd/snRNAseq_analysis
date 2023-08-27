#!/bin/bash
cd ~/snRNAseq_analysis/
mkdir output/patient_umaps/

for patient in {malignant,{A..P}} ; do
  echo $patient
  for umap in {infercnv,celltype} ; do
    cp \
      output/snRNAseq_workflow/by_patient_wo_organoids/$patient/integrating_mnn/integrating_mnn_files/figure-html/${umap}_annots_void_plot-1.png \
      output/patient_umaps/${patient}_${umap}_umap.png
  done
done



