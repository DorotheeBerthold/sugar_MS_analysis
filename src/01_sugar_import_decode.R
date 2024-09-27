###########################################################
## Metabolomic analysis of sugars in =MM12 timecourse.       
## Dorothée L. Berthold, ETH Zürich                
##########################################################

# Import libraries
{
  library(readxl)
  library(ggplot2)
  library(tidyverse)
  library(ggfortify)
  library(reshape2)
  library(remotes)
  library(DEGreport)
  library(pheatmap)
  library(vsn)
  library(ggforce)
  library(ggpubr)
  library(emmeans)
  library(ggforce)
  library(ggrepel)
  library(RColorBrewer)
  library(dendextend)
  library(openxlsx)
}
# Manage injection order & replicates
######################################################################################################################

skyline <- read.csv("tables/Andi_transitionlist.csv")
colnames(skyline)[1] <- "Molecule"

skyline$inject <- rownames(skyline)

# Find the first Peptide in the skyline data
first_peptide <- unique(skyline$Molecule)[1]

# Get the indices where the Peptide matches the first Peptide
inject_order <- which(skyline$Molecule == first_peptide)

# Repeat the injection order across all unique peptides
skyline$inject <- rep(inject_order, length(unique(skyline$Molecule)))


#Clean up table
skyline <- skyline %>% 
  select(Molecule, Area, Replicate, inject)

# Remove parentheses after standards
skyline$Replicate <- gsub("\\s*\\(.*\\)", "", skyline$Replicate)

# Concatenate Replicate.Name and inject columns with an underscore for pivoting
skyline$Replicate <- paste(skyline$Replicate, skyline$inject, sep = "_")


# Pivot wide
skyline_wide <- skyline%>% 
  pivot_wider(names_from = Molecule, values_from = Area, id_cols = Replicate)

# Convert all columns except the first to numeric and handle '#N/A' values
skyline_wide <- skyline_wide %>%
  mutate(across(-1, ~ as.numeric(na_if(., "#N/A")))) %>%
  # Split 'Replicate' column into 'Replicate' and 'inject'
  separate(Replicate, into = c("Replicate", "inject"), sep = "_", convert = TRUE)

# Import keys for sample decoding
keys <- read.csv("tables/conversion384-96well_2.csv")
colnames(keys) <- keys[1,]
colnames(keys)[1] <- "Replicate"
keys <-keys[-1,]

skyline_keys <- left_join(keys, skyline_wide, by = "Replicate")

# Import well to sample list for decoding
######################################################################################################################
code <- read.csv("tables/DB044_worklist_MS.csv")
code$string <- sub(".*_S\\d+\\s+(.+)$", "\\1", code$Data.File)

#filter out water & gc
code <- code %>% 
  filter(Sample.Type == "Sample")

code <- code[,c(3,9)]

meta <- read.csv("tables/DB044_meta.csv")
colnames(meta) [5] <- "string"

code <- code %>%
  separate(string, into = c("day", "hour", "replicate"), sep = "_", remove = F) 

code <- code %>% 
  mutate(hour = case_when(
    hour == "48h" ~ "24h",
    hour == "72h" ~ "24h",
    hour == "96h" ~ "24h",
    hour == "120h" ~ "24h",
    hour == "144h" ~ "24h",
    hour == "168h" ~ "24h",
    TRUE ~ hour),
    day = case_when(
      day == "media" ~ "blank",
      TRUE ~ day
    )
  )

code$string <- paste0(code$day, "_", code$hour, "_", code$replicate)
code <- code[,-c(3:5)]

meta_code <- left_join(code, meta, by = "string")

#Create new column based on position for matching with keys
meta_code$sample <- paste0("sample_", sub(".*-", "", meta_code$Sample.Position))


skyline_keys_meta <- full_join(meta_code, skyline_keys, by = "sample")
skyline_keys_meta <- skyline_keys_meta %>% 
  distinct(inject, .keep_all = TRUE) %>% 
  filter(!grepl("^H(0[5-9]|[1-9][0-9])$", `96well`))
skyline_keys_meta$sample <- gsub("^s([1-8])$", "std\\1", skyline_keys_meta$sample)


skyline_keys_meta <- skyline_keys_meta[,c(2,9:10,12:33)]

write.csv(skyline_keys_meta, "results/skyline_keys_meta.csv")
