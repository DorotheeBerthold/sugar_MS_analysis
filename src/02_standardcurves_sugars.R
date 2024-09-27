##########################################################
## Metabolomic analysis of sugars in OMM12 timecourse.   
## Dorothée L. Berthold, ETH Zürich                
#########################################################

# Manage injection order & replicates
######################################################################################################################

#standards <- read.csv("results/sugars_standards_decoded.csv", row.names = 1)
standards <- read.csv("results/skyline_keys_meta.csv", row.names = 1)


# Make table numeric (except first three columns)
#standards[,-(1:3) ] <- lapply(standards[,-(1:3) ], as.numeric)
standards[,-(1:4) ] <- lapply(standards[,-(1:4) ], as.numeric)

# Filter out mq

standards <- standards %>% 
  filter(sample != "mq")

#Define vector of all compounds
analytes <- c("Ribose", "Arabinose", "Xylose", "36anhydro", "Rhamnose",
              "Fucose", "Galactosamine", "Glucosamine", "Mannose", "Glucose",
              "Galactose", "13C-Glc", "13C-Man", "13C-Gal", "ManA", "GulA",
              "GlucuronicAcid", "GalacturonicAcid", "N-acetyl-Glucosamin",
              "N-acetyl-Mannosamin", "N-acetyl-Galactosamin")

current_names <- colnames(standards)

# Columns that need to be replaced
replace_names <- c("X36anhydro" = "36anhydro",
                   "X13C.Glc" = "13C-Glc",
                   "X13C.Man" = "13C-Man",
                   "X13C.Gal" = "13C-Gal",
                   "N.acetyl.Glucosamin" = "N-acetyl-Glucosamin",
                   "N.acetyl.Mannosamin" = "N-acetyl-Mannosamin",
                   "N.acetyl.Galactosamin" = "N-acetyl-Galactosamin")

# Replace column names only if they are in the replacement mapping
colnames(standards) <- ifelse(current_names %in% names(replace_names), 
                    replace_names[current_names], 
                    current_names)


# Normalize by calculating the mean of the three 13C-columns

standards <- standards %>%
  mutate(norm_inject = rowMeans(select(., `13C-Gal`, `13C-Man`, `13C-Glc`), na.rm = TRUE))


# Divide each sugar by norm_inject to see how stable the injections are

#standards_norm <- standards %>%
#  mutate(across(5:(ncol(standards)-1), ~ . / norm_inject))

# Normalize standards
######################################################################################################################
df_std <- standards%>% 
  filter(grepl("std", sample)) %>% 
  filter(!Replicate %in% c('P1-H03','P1-J01','P1-L01','P1-L05')) # These wells are apparently not great

mean_std <- mean(df_std$norm_inject)

df_std <- df_std %>% 
  mutate(norm_inject2 = norm_inject/mean_std)

df_std <- df_std %>%
  mutate(across(5:(ncol(df_std)-2), ~ . / norm_inject2))


# QC plot
######################################################################################################################

# Here, it appears highest standard is used as QC

qc <- df_std %>% 
  filter(grepl("std1", sample)) %>% 
  arrange(inject)

# Reshape data to long format
qc_long <- qc %>%
  select(inject, all_of(analytes)) %>%
  pivot_longer(cols = -inject, names_to = "analyte", values_to = "value")

# Create ggplot object with reshaped data
ggplot(qc_long, aes(x = inject, y = value, color = analyte)) +
  geom_line() +
  geom_point() +
  labs(
    title = "QC Analytes",
    x = "Injection Order",
    y = "Ion Counts"
  ) +
  theme(legend.position = "right")

ggsave("plots/QC_sugars.png", width = 8, height = 6, units = "in")

# Create standard curves
######################################################################################################################


# Extract level of standard dilution
df_std$level <- as.integer(substr(df_std$sample, 4, nchar(df_std$sample)))

# Calculate the amount based on the level:
df_std$amount <- as.numeric(3.33^(1 - df_std$level))

# Read in standard concentrations
std_conc <- read.csv("tables/Standards_sugars.csv")

# Create named list for standards (std_dict)
analytes <- std_conc$PythonName
std_dict <- setNames(as.list(std_conc$uM), analytes)

df_std <- df_std[,-c(1:3)]
df_std_long <- reshape2::melt(df_std, id.vars = c("level", "amount", "inject"))
colnames(df_std_long) [4:5] <- c("Molecule", "ratio")

# Filter out Molecules not found in std_dict
df_std_long <- df_std_long %>% 
  filter(Molecule %in% analytes)

# Multiply the amount by the concentration of each standard found in std_dict
df_std_long$concentration <- df_std_long$amount * 
  as.numeric(sapply(as.character(df_std_long$Molecule), function(x) std_dict[[x]]))

# Create linear model for slope fitting for each analyte
lm_std <- df_std_long %>% 
  group_by(Molecule) %>% 
  summarize(Slope = coef(lm(ratio ~ concentration))[2],
            Intercept = coef(lm(ratio ~ concentration))[1],
            R_squared = summary(lm(ratio~ concentration))$r.squared)

# Create a function to generate the equation text
lm_std <- lm_std %>%
  mutate(eq_label = paste0("y = ", round(Slope, 3), " * x + ", round(Intercept, 3), 
                           "\nR² = ", round(R_squared, 3)))




ggplot(df_std_long, aes(x = concentration, y = ratio, color = Molecule)) +
  geom_point() +
  geom_abline(data = lm_std, aes(intercept = Intercept, slope = Slope)) +
  facet_wrap(~Molecule, scales = "free_y") +
  theme_classic() +
  labs(title = "Standard curves sugars", x = "Concentration") +
  geom_text(data = lm_std, aes(x = Inf, y = Inf, label = eq_label, color = Molecule), 
            hjust = 1.1, vjust = 1.1, inherit.aes = F, size = 3)

ggsave("plots/standard_curves_sugars.png")

# Calculate back the concentrations from the ratios & the curves
######################################################################################################################

#Normalize samples
df_sample <- standards%>% 
  filter(grepl("sample", sample)) 

mean_samples <- mean(df_sample$norm_inject)

df_sample <- df_sample %>% 
  mutate(norm_inject2 = norm_inject/mean_std)

df_sample <- df_sample %>%
  mutate(across(5:(ncol(df_std)-2), ~ . / norm_inject2))

df_sample <- df_sample %>% 
  select(-c(norm_inject, norm_inject2, Replicate, inject))

# Make long df
df_sample_long <- reshape2::melt(df_sample, id.vars = c("string", "sample"))


# Split up string column
#df_sample_long <- df_sample_long %>% 
 # separate(string, into = c("day", "hour", "sample"), sep = "_", convert = TRUE)
colnames(df_sample_long) [3:4] <- c("Molecule", "ratio")

df_sample_long <- left_join(df_sample_long, lm_std, by = "Molecule")

df_sample_long <- df_sample_long %>%
  mutate(concentration = (ratio - Intercept) / Slope)

df_sample_long <- df_sample_long[,-c(4:8)]

# Pivot the data frame from long to wide format
df_wide <- df_sample_long %>%
  pivot_wider(
    id_cols = c(string, sample),    # Columns to keep as identifiers
    names_from = Molecule,             # Column to use for new column names
    values_from = concentration        # Column to use for values
  ) %>% 
  select(-c(`13C-Gal`, `13C-Man`, `13C-Glc`)) 

# Correct for dilution in sample prep (5uL in 100uL acid hydrolyis) & 27.5uL in NaOH (1.28)
df_wide[,-(1:2) ] <- lapply(df_wide[,-(1:2) ],function(x) x * 20*1.28) 
  

write.csv(df_wide, "results/sugars_conc_wide.csv")

