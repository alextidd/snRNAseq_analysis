#!/bin/bash
#PBS -N move
#PBS -o /rds/general/user/art4017/home/snRNAseq_analysis/log/move.out
#PBS -e /rds/general/user/art4017/home/snRNAseq_analysis/log/move.err
#PBS -l select=1:ncpus=1:mem=60gb
#PBS -l walltime=12:00:00

# tmux new -s snRNAseq_workflow

cd ~/snRNAseq_analysis/

mv output/snRNAseq_workflow/by_patient ~/../ephemeral/snRNAseq_workflow/output/