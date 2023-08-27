#!/bin/bash
counts=~/cellranger/output/
cd ~/snRNAseq_analysis
mkdir output
stats_file=~/snRNAseq_analysis/output/stats.tsv

echo -e "sample\tn_nuclei" > $stats_file ; 
(
  cd $counts
  cat sample_metadata.tsv |
  sed 1d |
  cut -f1 |
  sort -u |
  (
    while read sample ; do
      echo $sample
      n_nuclei=$(zcat $sample/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l)
      echo -e "$sample\t$n_nuclei" >> $stats_file
    done
  )
)