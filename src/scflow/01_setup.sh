#!/bin/bash

cd ~/snRNAseq_analysis/
data_dir=data/scflow/
sm_file=~/cellranger/output/sample_metadata.tsv

# create sample sheet
cat $sm_file |
sed 's/sample\t/manifest\t/g' |
awk -F'\t' 'OFS=FS{X=$1; gsub(/_/,".",X); $1=X}1' \
> $data_dir/sample_sheet.tsv

# create manifest file
echo -e 'key\tfilepath' > $data_dir/manifest.tsv
cat $sm_file |
sed 1d | cut -f1 |
(
  while read patient ; do
    echo $patient
    # '_' in manifest causes error, convert to '.'
    manifest=${patient//_/.}
    echo $manifest
    echo -e "$manifest\t/rds/general/user/art4017/home/cellranger/output/$patient/outs/filtered_feature_bc_matrix/" \
    >> $data_dir/manifest.tsv
  done
)
