####USING GENERAL SHEET in unedited NMDS####
library(tidyverse)
library(ggVennDiagram)
install.packages("vegan")
library(vegan)
library(stringr)

Bact_corrected <- Unedited_Bacteria %>%
  mutate(P1_Blank_Mean = rowMeans(across(c(P_B__1)), na.rm = TRUE)) %>%
  mutate(across(
    -c(ScientificName, P_S_1_1, P_S_2_1, P_S_3_1),
    ~ pmax(. - P1_Blank_Mean, 0)
  )) %>%
  mutate(P2_Blank_Mean = rowMeans(across(c(P_B__2)), na.rm = TRUE)) %>%
  mutate(across(
    -c(ScientificName, P_S_1_2, P_S_2_2, P_S_3_2),
    ~ pmax(. - P2_Blank_Mean, 0)
  )) %>%
  mutate(B1_Blank_Mean = rowMeans(across(c(B_B__1)), na.rm = TRUE)) %>%
  mutate(across(
    -c(ScientificName, B_S_1_1, B_S_2_1, B_S_3_1),
    ~ pmax(. - B1_Blank_Mean, 0)
  )) %>%
  mutate(B2_Blank_Mean = rowMeans(across(c(B_B__2)), na.rm = TRUE)) %>%
  mutate(across(
    -c(ScientificName, B_S_1_2, B_S_2_2, B_S_3_2),
    ~ pmax(. - B2_Blank_Mean, 0)
  )) %>%
  mutate(Blank_Mean = rowMeans(across(c(N_1, N_2)), na.rm = TRUE)) %>%
  mutate(across(
    -c(ScientificName, N_1, N_2, Blank_Mean),
    ~ pmax(. - Blank_Mean, 0)
  )) %>%
  mutate(across(
    starts_with("O_"),
    ~ pmax(. - O_B__1, 0)
  )) %>%
  mutate(across(
    starts_with("S_"),
    ~ pmax(. - S_B__2, 0)
  )) %>%
  select(-c(
    P_B__1, P_B__2, P1_Blank_Mean, P2_Blank_Mean,
    B_B__1, B_B__2, B1_Blank_Mean, B2_Blank_Mean,
    O_B__1, S_B__2, Blank_Mean, N_1, N_2
  ))

#  Standardise each row (species) to 0–1
Bact_rel <- Bact_corrected %>%
  rowwise() %>%
  mutate(across(-ScientificName, ~ {
    vals <- c_across(-ScientificName)
    rng <- range(vals, na.rm = TRUE)
    if (diff(rng) == 0) 0 else (. - rng[1]) / diff(rng)
  })) %>%
  ungroup()

########
# Remove ScientificName and transpose
rel_matrix <- Bact_rel %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame()

# Transpose so samples = rows, species = columns
rel_matrix_t <- t(rel_matrix)

# Remove samples with no data (row sum = 0)
rel_matrix_t <- rel_matrix_t[rowSums(rel_matrix_t) > 0, ]
set.seed(123)
nmds <- metaMDS(rel_matrix_t, distance = "bray", k = 2, trymax = 100)
nmds$stress

# Extract NMDS site scores
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Site <- rownames(site_scores)

site_scores <- site_scores %>%
  mutate(
    Site = str_replace_all(Site, " ", ""),
    Site = str_replace_all(Site, "__", "_"),
    SiteGroup = str_extract(Site, "^[A-Z]_[A-Z]")
  ) %>%
  filter(!is.na(SiteGroup))  %>%
  mutate(
    SiteGroup = factor(
      SiteGroup,
      levels = c("O_M","O_L","O_S","P_M","P_S", "B_S","B_M", "S_M","S_L","S_S")
    )
  ) 

# Extract first and second parts of SiteGroup
site_scores$Site <- substr(site_scores$SiteGroup, 1, 1)
site_scores$SampleType <- substr(site_scores$SiteGroup, 3, 3)

# Recode abbreviations to full names
site_scores$Site <- recode(site_scores$Site,
                           "O" = "Olifants",
                           "P" = "Platboom",
                           "B" = "Bordjiesrif",
                           "S" = "Smitswinkel")

site_scores$SampleType <- recode(site_scores$SampleType,
                                 "M" = "Metaprobe",
                                 "S" = "Sterivex",
                                 "L" = "Wilderlab")
site_scores$Site <- factor(site_scores$Site, levels = c(
  "Olifants", "Platboom", "Bordjiesrif","Smitswinkel"
))

metadata <- site_scores %>%
  select(Site, SampleType) %>%
  mutate(SampleID = rownames(site_scores))

# Match the order of metadata to the rows in rel_matrix_t
metadata <- metadata %>%
  filter(SampleID %in% rownames(rel_matrix_t)) %>%
  column_to_rownames("SampleID")

metadata <- metadata %>%
  filter(SampleID %in% rownames(rel_matrix_t)) %>%
  as.data.frame() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("SampleID")
all(rownames(rel_matrix_t) == rownames(metadata))
colnames(metadata)

library(dplyr)

temp_summary <- temp_summary %>%
  mutate(Site = recode(Site,
                       "O" = "Olifants",
                       "Plat" = "Platboom",
                       "Borj" = "Borjiesrif",
                       "Smits" = "Smits"))



metadata <- metadata %>%
  left_join(temp_summary %>% group_by(Site) %>% summarise(Temp_mean = mean(Temp_mean, na.rm = TRUE)),
            by = "Site")

set.seed(123)
adonis2(rel_matrix_t ~ Site + SampleType + Temp_mean,
        data = metadata,
        permutations = 999,
        method = "bray")

adonis2(rel_matrix_t ~ Temp_mean, data = metadata, permutations = 999, method = "bray")
adonis2(rel_matrix_t ~ Site, data = metadata, permutations = 999, method = "bray")
adonis2(rel_matrix_t ~ SampleType, data = metadata, permutations = 999, method = "bray")

#interaction
# Test interaction
adonis2(rel_matrix_t ~ Site * SampleType*Temp_mean, data = metadata, permutations = 999, method = "bray")


bray_dist <- vegdist(rel_matrix_t, method = "bray")

### 1. Test dispersion by Site
disp_site <- betadisper(bray_dist, group = metadata$Site)
anova(disp_site)  # Tests if dispersions are significantly different
permutest(disp_site)  # Permutation test

### 2. Test dispersion by Sampling Method
disp_method <- betadisper(bray_dist, group = metadata$SampleType)
anova(disp_method)
permutest(disp_method)

###Test disp by temp 
disp_temp <- betadisper(bray_dist, group = metadata$Temp_mean)
anova(disp_temp)  # Tests if dispersions are significantly different
permutest(disp_temp)
