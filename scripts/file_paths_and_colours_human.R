#.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-with_modules-4.4.0-3.sif", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths(c("/tgen_labs/jfryer/kolney/R/x86_64-pc-linux-gnu-library/4.3", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

.libPaths()

#----------------- Libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DESeq2) 
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(reshape2)
library(edgeR)  
library(limma)  
#library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(data.table)
library(philentropy)
library(gplots)
library(variancePartition)
library(NatParksPalettes) # colors


#BiocManager::install("variancePartition", lib="/tgen_labs/jfryer/YOURACCOUNTHERE/R/x86_64-pc-linux-gnu-library/4.3")
#BiocManager::install("variancePartition")
#devtools::install_github("DiseaseNeuroGenomics/variancePartition")
#BiocManager::install("limma", force = TRUE) 

#----------------- Define variables
#tissue <- c("Kidney") # Kidney or Brain
typeOfCount <- c("ReadsPerGene.out.tab") 
#pathToRef <- c("/tgen_labs/jfryer/projects/references/mouse/ensembl_v7/")

#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

color_disease <- c("gray", "green4")

#----------------- Data
metadata <- read.delim("/tgen_labs/jfryer/kolney/sepsis_human/metadata.tsv", header = TRUE, sep = "\t")

# samples to exclude 
# MDS Kidney outliers - F_40, M_96, M_22

# Four samples were sequenced on a separate lane 
# "102", "104", "106", "108",

# Cause of death was kidney failure 
# "56" - 10_82
samples_to_remove <- metadata %>%
  filter(!(banner_case_id %in% c("15_46", "12_20", "06_37", "10_82", "13_13", "13_22", "13_24"))) 
#          "94", "56", "32", "18", "4", "10", "46", "34")))

# Outliers on the MDS kidney plot: "15_46", "12_20", "06_37"

# control 46 - large sample specific weights
metadata <- samples_to_remove %>% # samples_to_keep
  mutate(
    sex_chr = if_else(sex == "F", "XX", "XY")
  )


metadata$Sample_ID <- paste0(metadata$tissue, "_", 
                             metadata$group, "_", 
                             metadata$sex, "_", 
                             metadata$unique_id)
