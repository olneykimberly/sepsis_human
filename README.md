# sepsis_human

Human brain and kidney samples from postmortem individuals. 

metadata <- read.delim("/tgen_labs/jfryer/kolney/sepsis_human/metadata/metadata.tsv")
sampleReadGroup <- read.delim("/tgen_labs/jfryer/kolney/sepsis_human/scripts/00_preprocessing/sampleReadGroupInfo.txt", header = FALSE)

merged_data <- merge(
  x = sampleReadGroup,
  y = metadata,
  by.x = "V1", 
  by.y = "R1_filename", 
  all = FALSE 
)


library(dplyr)
data_processed <- merged_data |>
  select(tissue, group, sex, unique_id, batch, V1, V2)
  
  
write.table(data_processed, "/tgen_labs/jfryer/kolney/sepsis_human/metadata/metadata_merged.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

