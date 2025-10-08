#!/bin/bash
#SBATCH --job-name=open_tar_sepsis_human                                                              
#SBATCH --time=10:00:00                               
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH -o slurm.open_tar_sepsis_human.job.%j.out
#SBATCH -e slurm.open_tar_sepsis_human.job.%j.err

## salloc --mem=64G --cpus-per-task=8
cd /tgen_labs/jfryer/arc_restore/
tar -xvzf human_banner_sepsis_bulk_RNAseq_raw_data.tar.gz