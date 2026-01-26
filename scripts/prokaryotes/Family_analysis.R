library(tidyverse)
library(dplyr)
library(tidyr)


####OTU ABUNDANCES ETC####
###using full dataset##
Correct <- Unedited_Bacteria %>%
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
df<-Correct

species_counts_per_site <- colSums(df[, -1] > 0)
min_nonzero_read <- min(df[, -1][df[, -1] > 0])
max_read <- max(Unedited_Bacteria[, -1])

# Print results
species_counts_per_site
min_nonzero_read
max_read

min_nonzero_read <- min(Unedited_Bacteria[, -1][Unedited_Bacteria[, -1] > 0])

species_counts_df <- data.frame(
  Site = names(species_counts_per_site),
  SpeciesCount = species_counts_per_site
)

####USING FAMILY SHEET####
Bact_corrected <- Bacteria %>%
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
  
####Standardise each row (species) to 0–1####
Bact_rel <- Bact_corrected %>%
  rowwise() %>%
  mutate(across(-ScientificName, ~ {
    vals <- c_across(-ScientificName)
    rng <- range(vals, na.rm = TRUE)
    if (diff(rng) == 0) 0 else (. - rng[1]) / diff(rng)
  })) %>%
  ungroup()

Bact_long_rel <- Bact_rel %>%
  pivot_longer(
    cols = -ScientificName,
    names_to = "Site",
    values_to = "RelAbundance"
  ) %>%
  mutate(Site = str_replace_all(Site, " ", ""))

Bact_long_rel <- Bact_long_rel %>%
  mutate(Group = substr(Site, 1, 1))

site_order <- Bact_long_rel %>%
  distinct(Site) %>%
  mutate(
    Site_chr = as.character(Site),
    Prefix = substr(Site_chr, 1, 1),
    SampleType = str_extract(Site_chr, "_[A-Z]+_?"), 
    Suffix = substr(Site_chr, nchar(Site_chr), nchar(Site_chr)),
    SampleType = recode(SampleType, `_B__` = "_B_"),  
    SampleOrder = factor(SampleType, levels = c("_L_", "_S_", "_M_"), ordered = TRUE)
  ) %>%
  arrange(
    factor(Prefix, levels = c("O", "P", "B", "S", "N")),
    SampleOrder,
    Suffix
  ) %>%
  pull(Site)




Bact_long_rel$Site <- factor(Bact_long_rel$Site, levels = site_order)
Bact_long_rel$Group <- factor(Bact_long_rel$Group, levels = c("O", "P", "B", "S"))

ggplot(Bact_long_rel, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(
    x = "Site",
    y = "Relative Abundance"
  )

length(unique(Bact_long_rel$ScientificName))
####SPLIT ACTIVE####
####SEPARATE FOR ACTIVE AND PASSIVE####
Bact_long_rel <- Bact_long_rel %>%
  mutate(SamplingType = case_when(
    str_detect(Site, "_M_") ~ "Passive",
    str_detect(Site, "_L_") ~ "Wilderlab",
    str_detect(Site, "_S_") ~ "Sterivex",
    TRUE ~ "Other"
  ))

####WILDERLAB####
wilder_data <- Bact_long_rel %>%
  filter(SamplingType == "Wilderlab") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Bact_long_rel %>%
  filter(RelAbundance > 0) %>%
  group_by(Site) %>%
  summarise(FamilyCount = n_distinct(ScientificName))

wilder_family_counts <- family_counts %>%
  filter(Site %in% wilder_data$Site) %>%
  mutate(Site = factor(Site))

wilder_data <- wilder_data %>%
  left_join(wilder_family_counts, by = "Site")

species_counts_df_wilder <- species_counts_df %>%
  filter(Site %in% unique(wilder_data$Site)) %>%
  mutate(Site = factor(Site, levels = levels(wilder_data$Site)))
label_data <- left_join(wilder_family_counts, species_counts_df_wilder, by = "Site")

site_levels <- levels(wilder_data$Site)
site_levels
pos1 <- which(site_levels == "O_L_3_2")
pos2 <- which(site_levels == "S_L_1_1")
midpoint <- (pos1 + pos2) / 2

# Calculate total abundance per site for label positioning
total_abundance_per_site <- wilder_data %>%
  group_by(Site) %>%
  summarise(TotalAbundance = sum(RelAbundance))

# Plot
p_wilder <- ggplot(wilder_data, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
  geom_text(
    data = label_data %>% left_join(total_abundance_per_site, by = "Site"),
    aes(x = Site, y = 1.02 * TotalAbundance, label = FamilyCount),
    inherit.aes = FALSE,
    size = 7,
    vjust = 0
  ) +
  geom_text(
    data = label_data,
    aes(x = Site, y = 88, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 7,
    vjust = 0
  ) +
  
  geom_vline(xintercept = midpoint, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),  # X-axis title
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(
    x = "Sample",
    y = "Relative Read Abundance",
    fill = "Family"
  )
p_wilder

####STERIVEX####
steri_data <- Bact_long_rel %>%
  filter(SamplingType == "Sterivex") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Bact_long_rel %>%
  filter(RelAbundance > 0) %>%
  group_by(Site) %>%
  summarise(FamilyCount = n_distinct(ScientificName))

steri_family_counts <- family_counts %>%
  filter(Site %in% steri_data$Site) %>%
  mutate(Site = factor(Site))

steri_data <- steri_data %>%
  left_join(steri_family_counts, by = "Site")

species_counts_df_steri <- species_counts_df %>%
  filter(Site %in% unique(steri_data$Site)) %>%
  mutate(Site = factor(Site, levels = levels(steri_data$Site)))
label_data <- left_join(steri_family_counts, species_counts_df_steri, by = "Site")

site_levels <- levels(steri_data$Site)
site_levels
pos1 <- which(site_levels == "P_S_3_2")
pos2 <- which(site_levels == "B_S_1_1")
midpoint <- (pos1 + pos2) / 2

# Calculate total abundance per site for label positioning
total_abundance_per_site <- steri_data %>%
  group_by(Site) %>%
  summarise(TotalAbundance = sum(RelAbundance))

# Plot
p_steri <- ggplot(steri_data, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
  geom_text(
    data = label_data %>% left_join(total_abundance_per_site, by = "Site"),
    aes(x = Site, y = 1.02 * TotalAbundance, label = FamilyCount),
    inherit.aes = FALSE,
    size = 7,
    vjust = 0
  ) +
  geom_text(
    data = label_data,
    aes(x = Site, y = 78, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  
  geom_vline(xintercept = midpoint, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),  # X-axis title
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 9),
    legend.position = "none"
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(
    x = "Sample",
    y = "Relative Read Abundance",
    fill = "Family"
  )
p_steri

####SEPARATE FOR ACTIVE AND PASSIVE####
Bact_long_rel <- Bact_long_rel %>%
  mutate(SamplingType = case_when(
    str_detect(Site, "_M_") ~ "Passive",
    str_detect(Site, "(_L_|_S_)") ~ "Active",
    TRUE ~ "Other"
  ))
table(Bact_long_rel$SamplingType)

active_data <- Bact_long_rel %>%
  filter(SamplingType == "Active") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Bact_long_rel %>%
  filter(RelAbundance > 0) %>%
  group_by(Site) %>%
  summarise(FamilyCount = n_distinct(ScientificName))

active_family_counts <- family_counts %>%
  filter(Site %in% active_data$Site) %>%
  mutate(Site = factor(Site))

active_data <- active_data %>%
  left_join(active_family_counts, by = "Site")

species_counts_df_active <- species_counts_df %>%
  filter(Site %in% unique(active_data$Site)) %>%
  mutate(Site = factor(Site, levels = levels(active_data$Site)))
label_data <- left_join(active_family_counts, species_counts_df_active, by = "Site")

site_levels <- levels(active_data$Site)
pos1 <- which(site_levels == "P_S_3_2")
pos2 <- which(site_levels == "B_S_1_1")
midpoint <- (pos1 + pos2) / 2

# Calculate total abundance per site for label positioning
total_abundance_per_site <- active_data %>%
  group_by(Site) %>%
  summarise(TotalAbundance = sum(RelAbundance))

# Plot
p_active <- ggplot(active_data, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
  geom_text(
    data = label_data %>% left_join(total_abundance_per_site, by = "Site"),
    aes(x = Site, y = 1.02 * TotalAbundance, label = FamilyCount),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  geom_text(
    data = label_data,
    aes(x = Site, y = 85, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  
  geom_vline(xintercept = midpoint, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Site",
    y = "Relative Abundance"
  )+
  theme(legend.position = "none")
p_active


passive_data <- Bact_long_rel %>%
  filter(SamplingType == "Passive") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )


family_counts <- Bact_long_rel %>%
  filter(RelAbundance > 0) %>%
  group_by(Site) %>%
  summarise(FamilyCount = n_distinct(ScientificName))

passive_family_counts <- family_counts %>%
  filter(Site %in% passive_data$Site) %>%
  mutate(Site = factor(Site))

passive_data <- passive_data %>%
  left_join(passive_family_counts, by = "Site")

species_counts_df_passive <- species_counts_df %>%
  filter(Site %in% unique(passive_data$Site)) %>%
  mutate(Site = factor(Site, levels = levels(passive_data$Site)))

label_data <- left_join(passive_family_counts, species_counts_df_passive, by = "Site")

site_levels <- levels(passive_data$Site)
pos1 <- which(site_levels == "P_M_3_2")
pos2 <- which(site_levels == "B_M_1_2")
midpoint <- (pos1 + pos2) / 2

# Calculate total abundance per site for label positioning
total_abundance_per_site <- passive_data %>%
  group_by(Site) %>%
  summarise(TotalAbundance = sum(RelAbundance))

# Plot
p_passive <- ggplot(passive_data, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
  geom_text(
    data = label_data %>% left_join(total_abundance_per_site, by = "Site"),
    aes(x = Site, y = 1.02 * TotalAbundance, label = FamilyCount),
    inherit.aes = FALSE,
    size = 7,
    vjust = 0
  ) +
  geom_text(
    data = label_data,
    aes(x = Site, y = 90, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 7,
    vjust = 0
  ) +
  
  geom_vline(xintercept = midpoint, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),  # X-axis title
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(
    x = "Sample",
    y = "Relative Read Abundance")+
      theme(legend.position = "none")
p_passive


##legend 
library(cowplot)

# Extract just the legend
legend <- get_legend(
  p_wilder + theme(legend.position = "right")
)

# Show only the legend
plot(legend)



####USING PRIMERS####
Primers<-Bacteria

Primers[is.na(Primers)] <- 0
Primers[Primers == ""] <- 0

Primers[, -1] <- lapply(Primers[, -1], as.numeric)

print(Primers)

primers_16S <- Primers$Name[Primers$UM == 1]
primers_18S <- Primers$Name[Primers$BU == 1]

venn_data <- list('16S' = primers_16S, '18S' = primers_18S)
install.packages("ggVennDiagram")
library(ggVennDiagram)

####Gradient####
ggVennDiagram(venn_data, label = "count") +
  theme_void() +                         
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA) )
ggplot2::theme_minimal()

####Jaccard####
presence_matrix <- Primers[, -1]  # remove 'Species' column
presence_matrix[is.na(presence_matrix)] <- 0
#transpose
t_matrix <- t(presence_matrix)
library(vegan)
jaccard_dist <- vegdist(t_matrix, method = "jaccard", binary = TRUE)
jaccard_similarity <- 1 - as.matrix(jaccard_dist)
print(round(jaccard_similarity, 2))

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

#### Standardise each row to 0–1 ####
Bact_rel <- Bact_corrected %>%
  rowwise() %>%
  mutate(across(-ScientificName, ~ {
    vals <- c_across(-ScientificName)
    rng <- range(vals, na.rm = TRUE)
    if (diff(rng) == 0) 0 else (. - rng[1]) / diff(rng)
  })) %>%
  ungroup()

rel_matrix <- Bact_rel %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame()
rel_matrix_t <- t(rel_matrix)

# Remove samples with no data (row sum = 0)
rel_matrix_t <- rel_matrix_t[rowSums(rel_matrix_t) > 0, ]
set.seed(123)
nmds <- metaMDS(rel_matrix_t, distance = "bray", k = 2, trymax = 100)
# Extract NMDS site scores
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Site <- rownames(site_scores)

#Cleaning and extraction
site_scores <- site_scores %>%
  mutate(
    Site = str_replace_all(Site, " ", ""),
    Site = str_replace_all(Site, "__", "_"),
    SiteGroup = str_extract(Site, "^[A-Z]_[A-Z]")
  ) %>%
  filter(!is.na(SiteGroup))  

site_scores$Site <- substr(site_scores$SiteGroup, 1, 1)
site_scores$SampleType <- substr(site_scores$SiteGroup, 3, 3)

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

site_scores$NMDS1 <- -site_scores$NMDS1
site_scores$NMDS2 <- -site_scores$NMDS2
ggplot(site_scores, aes(x = NMDS1, y = -NMDS2,
                        color = Site,
                        shape = SampleType)) +
  geom_point(size = 6) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_color_manual(values = c(
    "Olifants" = "blue",
    "Platboom" = "orange",
    "Bordjiesrif" = "#CC79A7",
    "Smitswinkel" = "green"
  )) +
  scale_shape_manual(values = c(
    "Metaprobe" = 16,   # circle
    "Sterivex" = 17,    # triangle
    "Wilderlab" = 15    # square
  )) +
  labs(   x = "NMDS1",
          y = "NMDS2",
          color = "Site",
          shape = "Sample Type"
  )


# Match metadata to matrix
meta <- site_scores %>%
  select(Site, SiteGroup) %>%
  filter(Site %in% rownames(rel_matrix_t)) %>%
  arrange(match(Site, rownames(rel_matrix_t)))

Site <- rownames(rel_matrix_t)
Site <- Site %>%
  str_replace_all(" ", "") %>%
  str_replace_all("__", "_")
meta <- data.frame(
  Site = Site,
  SiteGroup = str_extract(Site, "^[A-Z]_[A-Z]")
) %>%
  filter(!is.na(SiteGroup))  
rownames(meta) <- meta$Site
rel_matrix_t <- rel_matrix_t[rownames(meta), ]  

# Run PERMANOVA
adonis_result <- adonis2(rel_matrix_t ~ SiteGroup, data = meta, method = "bray", permutations = 999)
print(adonis_result)
disp <- betadisper(vegdist(rel_matrix_t, method = "bray"), meta$SiteGroup)
anova(disp)  # significance test

####PASSIVE SAMPLES####
# Select only "_M" columns 
metaprobe_cols <- names(Bact_rel)[str_detect(names(Bact_rel), "_M")]

df_meta <- Bact_rel %>%
  select(ScientificName, all_of(metaprobe_cols))

meta_matrix <- df_meta %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame() %>%
  t()  

# Remove samples with all-zero rows
meta_matrix <- meta_matrix[rowSums(meta_matrix) > 0, ]

set.seed(123)
nmds_meta <- metaMDS(meta_matrix, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores
meta_scores <- as.data.frame(scores(nmds_meta, display = "sites"))
meta_scores$Site <- rownames(meta_scores)

# Clean site names and extract SiteGroup 
meta_scores <- meta_scores %>%
  mutate(
    Site = str_replace_all(Site, " ", ""),
    Site = str_replace_all(Site, "__", "_"),
    SiteGroup = str_extract(Site, "^[A-Z]_[A-Z]")
  ) %>%
  filter(!is.na(SiteGroup))  %>%
  mutate(
    SiteGroup = factor(
      SiteGroup,
      levels = c("O_M","P_M", "B_M", "S_M")
    )
  ) 

ggplot(meta_scores, aes(x = NMDS1, y = NMDS2, color = SiteGroup, label = Site)) +
  geom_point(size = 4) +
  geom_text(vjust = -1.2, size = 3) +
  theme_minimal() +
  scale_color_discrete(
    labels = c("O_M" = "Olifants Metaprobe", "P_M" = "Platboom Metaprobe",
               "B_M" = "Bordjiesrif Metaprobe", "S_M" = "Smitswinkel Metaprobe"))+
  labs(
    title = "NMDS of Metaprobe Samples (Relative Abundance)",
    x = "NMDS1",
    y = "NMDS2",
    color = "Site Group"
  )

# Prepare metadata
meta <- meta_scores %>%
  select(Site, SiteGroup) %>%
  filter(Site %in% rownames(meta_matrix)) %>%
  arrange(match(Site, rownames(meta_matrix)))

stopifnot(all(meta$Site == rownames(meta_matrix)))

# Run PERMANOVA
adonis_meta <- adonis2(meta_matrix ~ SiteGroup, data = meta, method = "bray", permutations = 999)
print(adonis_meta)
disp <- betadisper(vegdist(meta_matrix, method = "bray"), meta$SiteGroup)
anova(disp) 

####ACTIVE####
active_cols <- names(Bact_rel)[str_detect(names(Bact_rel), "_S|_L")]

df_active <- Bact_rel %>%
  select(ScientificName, all_of(active_cols))
active_matrix <- df_active %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame() %>%
  t()  

# Remove samples with all-zero rows
active_matrix <- active_matrix[rowSums(active_matrix) > 0, ]

set.seed(123)
nmds_active <- metaMDS(active_matrix, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores
active_scores <- as.data.frame(scores(nmds_active, display = "sites"))
active_scores$Site <- rownames(active_scores)

# Clean site names and extract SiteGroup
active_scores <- active_scores %>%
  mutate(
    Site = str_replace_all(Site, " ", ""),
    Site = str_replace_all(Site, "__", "_"),
    SiteGroup = str_extract(Site, "^[A-Z]_[A-Z]")
  ) %>%
  filter(!is.na(SiteGroup)) %>%
  mutate(
    SiteGroup = factor(
      SiteGroup,
      levels = c("O_S","O_L", "P_S", "B_S", "S_S","S_L")
    )
  ) 

ggplot(active_scores, aes(x = NMDS1, y = NMDS2, color = SiteGroup, label = Site)) +
  geom_point(size = 4) +
  geom_text(vjust = -1.2, size = 3) +
  theme_minimal() +
  scale_color_discrete(
    labels = c("O_S" = "Olifants Sterivex", "P_S" = "Platboom Sterivex",
               "B_S" = "Bordjiesrif Sterivex", "S_S" = "Smitswinkel Sterivex", "S_L" = "Smitswinkel Wilderlab","O_L" = "Olifants Wilderlab"))+
  labs(
    title = "NMDS of Sterivex Samples (Relative Abundance)",
    x = "NMDS1",
    y = "NMDS2",
    color = "Site Group"
  )

# Prepare metadata
active <- active_scores %>%
  select(Site, SiteGroup) %>%
  filter(Site %in% rownames(active_matrix)) %>%
  arrange(match(Site, rownames(active_matrix)))

stopifnot(all(active$Site == rownames(active_matrix)))

# Run PERMANOVA
adonis_active <- adonis2(active_matrix ~ SiteGroup, data = active, method = "bray", permutations = 999)
print(adonis_active)
disp <- betadisper(vegdist(active_matrix, method = "bray"), active$SiteGroup)
anova(disp)


####Indicator Species####
sample_names <- rownames(rel_matrix_t)

# Define site groups
site_group <- case_when(
  str_detect(sample_names, "O_") ~ "Olifants",
  str_detect(sample_names, "P_") ~ "Platboom",
  str_detect(sample_names, "B_") ~ "Bordjiesrif",
  str_detect(sample_names, "S_") ~ "Smitswinkel",
  TRUE ~ "Other"
)

library(indicspecies)

# Run indicator species analysis
set.seed(123)
isa_sites <- multipatt(rel_matrix_t, site_group, func = "IndVal.g", control = how(nperm = 999))

# Summary
summary(isa_sites, indvalcomp = TRUE, alpha = 0.05)

group_labels <- c("Bordjiesrif", "Olifants", "Platboom", "Smitswinkel")

library(dplyr)
library(tibble)
sign_df <- as.data.frame(isa_sites$sign) %>% 
  rownames_to_column("Species")
stat_df <- as.data.frame(isa_sites$stat) %>%
  rownames_to_column("Species")

# Merge them by species
results_df <- left_join(sign_df, stat_df, by = "Species")

library(dplyr)
library(gt)

# Filter significant indicators (p.value < 0.05)
significant_results <- results_df %>%
  filter(p.value < 0.05)

# Identify which site each species is associated with
# If 1 indicates presence and 0 absence
significant_results <- significant_results %>%
  mutate(Associated_Sites = paste0(
    ifelse(s.Bordjiesrif == 1, "Bordjiesrif, ", ""),
    ifelse(s.Olifants == 1, "Olifants, ", ""),
    ifelse(s.Platboom == 1, "Platboom, ", ""),
    ifelse(s.Smitswinkel == 1, "Smitswinkel, ", "")
  )) %>%
  mutate(Associated_Sites = gsub(", $", "", Associated_Sites)) 

summary_table <- significant_results %>%
  select(Species, Associated_Sites, stat, p.value) %>%
  arrange(Associated_Sites, p.value)

my_gt_table <-summary_table %>%
  gt() %>%
  tab_header(
    title = "Significant Indicator Species by Site"
  ) %>%
  fmt_number(
    columns = c(stat, p.value),
    decimals = 3
  ) %>%
  cols_label(
    Species = "Taxon",
    Associated_Sites = "Indicator Site(s)",
    stat = "Indicator Value",
    p.value = "p-value"
  ) %>%
  data_color(
    columns = vars(p.value),
    colors = scales::col_numeric(
      palette = c("green", "yellow", "red"),
      domain = c(0, 0.05)
    )
  )

gtsave(my_gt_table, "bact_indicator_species_table.png")


####retry with black labels####
#### REWORKED INDICATOR SPECIES TREE ####
#### FULL REWORKED INDICATOR SPECIES TREE ####
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(ggplot2)

indicator_species <- summary_table$Species   

indicators_only <- rel_matrix_t[, colnames(rel_matrix_t) %in% indicator_species]

# Remove species with all zero counts
indicators_only <- indicators_only[, colSums(indicators_only) > 0]
# Build the OTU clustering tree 
otu_dist <- vegdist(t(indicators_only), method = "bray")
otu_hclust <- hclust(otu_dist, method = "average")
otu_tree <- as.phylo(otu_hclust)

tip_metadata <- data.frame(
  label = otu_tree$tip.label,
  Associated_Sites = significant_results$Associated_Sites[
    match(otu_tree$tip.label, significant_results$Species)
  ],
  stringsAsFactors = FALSE
)

tip_metadata$Associated_Sites[is.na(tip_metadata$Associated_Sites)] <- "None"

tip_metadata_long <- tip_metadata %>%
  separate_rows(Associated_Sites, sep = ",\\s*")  

# Map tip labels to numeric y positions and create stack offset
tip_metadata_long <- tip_metadata_long %>%
  left_join(data.frame(label = otu_tree$tip.label, y = 1:length(otu_tree$tip.label)), by = "label") %>%
  group_by(label) %>%
  mutate(circle_offset = seq_along(Associated_Sites) - 1)  

# Set custom legend order 
all_categories <- unique(tip_metadata_long$Associated_Sites)
priority_cats <- c("Olifants", "Platboom","Olifants, Platboom","Bordjiesrif", "Smitswinkel", "Bordjiesrif, Smitswinkel")
other_cats <- setdiff(all_categories, priority_cats)
final_levels <- c(priority_cats, sort(other_cats))
tip_metadata_long$Associated_Sites <- factor(tip_metadata_long$Associated_Sites, levels = final_levels)

# Assign colors
manual_colors <- c(
  "Bordjiesrif" = "#CC79A7",
  "Olifants" = "blue",
  "Platboom" = "orange",
  "Smitswinkel" = "green",
  "Olifants, Platboom"= "yellow",
  "Bordjiesrif, Smitswinkel"="red"
)

auto_cats <- setdiff(final_levels, names(manual_colors))

auto_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(auto_cats)),
  auto_cats
)

combined_colors <- c(manual_colors, auto_colors)

# Ensure 'None' is black
if("None" %in% names(combined_colors)) combined_colors["None"] <- "black"


library(dplyr)
library(tidyr)
library(ggtree)

# Split species with multiple sites into separate rows
tip_metadata_long <- tip_metadata %>%
  separate_rows(Associated_Sites, sep = ",\\s*")  

# Add stacking for multiple sites
tip_metadata_long <- tip_metadata_long %>%
  group_by(label) %>%
  mutate(stack = row_number() - 1) %>% 
  ungroup()

# Plot circular tree
tree_plot <- ggtree(otu_tree, layout = "circular") %<+% tip_metadata_long +
  geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.02, linesize = 0.2) +
  geom_point2(aes(color = Associated_Sites),
              size = 2.5,
              offset = 0.02,
              nudge = 0.02 * tip_metadata_long$stack) +  # stack circles for multi-site
  scale_color_manual(values = combined_colors, name = "Indicator Site(s)") +
  theme_void() +
  theme(legend.position = "right")

tree_plot

###rearrange tips
library(ggtree)
library(dplyr)
library(tidyr)

# Split multi-site species into separate rows
tip_metadata_long <- tip_metadata %>%
  separate_rows(Associated_Sites, sep = ",\\s*") %>%
  group_by(label) %>%
  mutate(stack = row_number() - 1) %>%
  ungroup()

# Base tree
tree_base <- ggtree(otu_tree, layout = "circular")

#Extract tip coordinates
tip_coords <- tree_base$data %>%
  filter(isTip) %>%
  select(label, x, y)

# Join metadata with coordinates 
plot_data <- tip_metadata_long %>%
  left_join(tip_coords, by = "label") %>%
  mutate(x = x + 0.02 + 0.02 * stack) 

# Plot 
tree_plot <- tree_base +
  geom_tiplab(aes(label = label),
              size = 2,
              align = TRUE,
              offset = 0.08,   
              linesize = 0.2) +
  geom_point(data = plot_data,
             aes(x = x, y = y, color = Associated_Sites),
             size = 2) +
  scale_color_manual(values = combined_colors, name = "Indicator Site(s)") +
  theme_void() +
  theme(legend.position = "right")

tree_plot
ggsave(
  filename = "tree_plot.png",
  plot = tree_plot,
  width = 16,       
  height = 17,      
  units = "in",
  dpi = 1200,
  limitsize = FALSE
)


