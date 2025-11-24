# sepsis_human
The transcriptional alternations within the kidney associated with sepsis.

| Group                     | Count   |
| ------------------------- |:-------:|
|human kidney controls      | 17      | 
|human kidney sepsis        | 9       | 
|pig kidney E.coli          | 4      | 
|pig kidney saline          | 4       | 

This git repo contains scripts for the following:
-   Metadata analysis
-   Processing of bulk RNA-sequencing data
-   Generation of manuscript figures 


## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the conda environment:
```
conda env create -n sepsis_human --file sepsis_human.yml

# To activate this environment, use
#
#     $ conda activate sepsis_human
#
# To deactivate an active environment, use
#
#     $ conda deactivate sepsis_human
```

## Reference genome
Reference genome and annotation were downloaded prior to running snakemake. 
1. human - Gencode GRCh38
2. pig - Ensembl Sus_scrofa.Sscrofa11.1
```
# 1. human
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

# 2. pig
cd reference
wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-107/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.103.gtf.gz 
```

## Snakemake for trimming reads and alignment to the reference genome
See the snakefile for specific details for each step. 
The config file was generated using the 01 and 02 scripts for obtaining sample information and creating the config file. 
```
Snakemake -s Snakefile 
```
Output includes STAR reads per gene, which can then be read into R for differential expression.

## Variance assessment and differential expression
```
R 01_human_differential_expression.Rmd
R 02_pig_differential_expression.Rmd
```
Output is differentially expressed genes for human sepsis vs control and for pig E.coli vs saline comparisons.

## References
All packages used in this workflow are publicly available. If you use this workflow please cite the packages used. 
If you use the data in this workflow cite the following:
[et al. 2026]()

## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | kolney@tgen.org |
| John Fryer, PhD | jfryer@tgen.org |

