#!/bin/bash

cd ~/snRNAseq_analysis/
data_dir=data/scflow/
sm_file=~/cellranger/output/sample_metadata.tsv

# create sample manifest file
(
  echo -e 'id\tid_col\tdir' ;
  cat $sm_file |
  awk -F'\t' 'NR>1{print $1"\tsample\t/rds/general/user/art4017/home/cellranger/output/"$1"/outs/filtered_feature_bc_matrix/"}'
) | cat > data/snRNAseq_workflow/sample_manifest.tsv