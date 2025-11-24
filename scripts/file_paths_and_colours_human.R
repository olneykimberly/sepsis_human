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
library(reshape2)
library(dplyr)
library(biomaRt)
library(WGCNA)
library(gprofiler2)
library(ComplexUpset)
library(patchwork)


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

plotCorrMatrix_gg <- function(C,
                              title = "Canonical Correlation Matrix",
                              reorder = FALSE,   # TRUE to reorder by hierarchical clustering
                              show_values = TRUE,
                              value_digits = 2) {
  
  # Input validation / conversion ------------------------------------------------
  if (is.null(dim(C))) stop("C must be a matrix or data.frame-like object.")
  # Force to matrix
  Cmat <- as.matrix(C)
  
  # Try to coerce to numeric, warn if any non-numeric entries
  if (!is.numeric(Cmat)) {
    warning("Coercing correlation matrix values to numeric.")
    Cmat <- apply(Cmat, c(1,2), function(x) as.numeric(as.character(x)))
  }
  
  if (any(is.na(Cmat))) {
    warning("NAs present in the correlation matrix after coercion.")
  }
  
  # Optional reordering by hierarchical clustering -------------------------------
  if (isTRUE(reorder)) {
    # distance based on 1 - abs(corr) so correlated items cluster together
    d <- as.dist(1 - abs(Cmat))
    hc <- hclust(d)
    ord <- hc$labels[hc$order]
    Cmat <- Cmat[ord, ord]
  }
  
  # Convert to long format (safe; avoids reshape2/data.table melt issues) -------
  df <- as.data.frame(as.table(Cmat))
  names(df) <- c("Var1", "Var2", "Correlation")
  
  # Keep factor levels in the matrix order
  df$Var1 <- factor(df$Var1, levels = rev(rownames(Cmat)))
  df$Var2 <- factor(df$Var2, levels = colnames(Cmat))
  
  # Build ggplot ----------------------------------------------------------------
  p <- ggplot(df, aes(x = Var2, y = Var1, fill = Correlation)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(min(df$Correlation, na.rm = TRUE),
                               max(df$Correlation, na.rm = TRUE)),
      name = "Corr"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )
  
  if (isTRUE(show_values)) {
    p <- p + geom_text(aes(label = sprintf(paste0("%.", value_digits, "f"), Correlation)),
                       size = 2.8)
  }
  
  return(p)
}


fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

makePaddedDataFrame <- function(l, ...) {
  maxlen <- max(sapply(l, length))
  data.frame(lapply(l, na.pad, len = maxlen), ...)
}


color_disease <- c("gray50", "brown4")

color_disease_human <- c("gray30", "brown4")
color_disease_pig <- c("gray30", "brown4")

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
