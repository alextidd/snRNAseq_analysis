#!/bin/bash
data_dir=~/cellranger/output/
sample_metadata=~/cellranger/output/sample_metadata.tsv
wd=~/snRNAseq_analysis/
script_dir=$wd/src/qsub/
vhl_dir=~/vhl_snRNAseq/
mkdir -p $script_dir
cd ~/snRNAseq_analysis/

cores=2
mem=30

# run generate_qc_report

cat $sample_metadata |
sed 1d |
cut -f1 |
sort -u |
(
  while read sample ; do
    echo $sample
    run_info=${sample}
    script=$script_dir/${run_info}.qsub
    
    > ${wd}/log/qsub_${run_info}.err
    ${wd}/log/qsub_${run_info}.out
    
cat << EOF > $script
#PBS -N ${run_info}

#PBS -l select=1:ncpus=${cores}:mem=${mem}gb
#PBS -l walltime=24:00:00
#PBS -o ${wd}/log/qsub_${run_info}.out
#PBS -e ${wd}/log/qsub_${run_info}.err

echo "running vhl_snRNAseq ${sample} from the oesophageal 10X data"

# codify params
sample=${sample}

# get directories
cd $wd

# conda env
. ~/.bashrc
conda activate vhl2

# run generate_QC_report.R
Rscript $wd/src/run_vhl.R \
  $data_dir/$sample/outs/filtered_feature_bc_matrix/ \
  $wd/output/$sample/ \
  $sample
EOF

      # submit the script
      qsub $script
      
    done
)