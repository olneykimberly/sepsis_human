#!/bin/bash
#SBATCH --job-name=human_sepsis                                                             
#SBATCH --time=36:00:00                               
#SBATCH --mem=1G
#SBATCH -n 1 # threaded 
#SBATCH -o slurm.human_sepsis.job.%j.out
#SBATCH -e slurm.human_sepsis.job.%j.err

# activate conda environment
source $HOME/.bash_profile
#module load python3
conda activate sepsis_human

# change directory to where Snakefile is located
CWD="/tgen_labs/jfryer/kolney/sepsis_human/scripts/00_preprocessing"
cd $CWD

snakemake --nolock -s Snakefile --jobs 86 --executor slurm --profile slurm_profile --rerun-incomplete --default-resources mem_mb=64000 ntasks=1 runtime=60 cpus_per_task=8