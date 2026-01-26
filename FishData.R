library(tidyverse)
library(dplyr)
library(tidyr)

####OTU ABUNDANCES ETC####
###using full dataset##
Full <- Full_Fish_Dataset %>%
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
df<-Full

species_counts_per_site <- colSums(df[, -1] > 0)
min_nonzero_read <- min(df[, -1][df[, -1] > 0])
max_read <- max(df[, -1])

# Print results
species_counts_per_site
min_nonzero_read
max_read

species_counts_df <- data.frame(
  Site = names(species_counts_per_site),
  SpeciesCount = species_counts_per_site
)
####USING FAMILY SHEET in edited dataset####
Fish_corrected <- Fish %>%
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

# Standardise each row (species) to 0–1
Fish_rel <- Fish_corrected %>%
  rowwise() %>%
  mutate(across(-ScientificName, ~ {
    vals <- c_across(-ScientificName)
    rng <- range(vals, na.rm = TRUE)
    if (diff(rng) == 0) 0 else (. - rng[1]) / diff(rng)
  })) %>%
  ungroup()

Fish_long_rel <- Fish_rel %>%
  pivot_longer(
    cols = -ScientificName,
    names_to = "Site",
    values_to = "RelAbundance"
  ) %>%
  mutate(Site = str_replace_all(Site, " ", ""))

Fish_long_rel <- Fish_long_rel %>%
  mutate(Group = substr(Site, 1, 1))

# Custom site ordering
site_order <- Fish_long_rel %>%
  distinct(Site) %>%
  arrange(factor(substr(Site, 1, 1), levels = c("O", "P", "B", "S", "N"))) %>%
  pull(Site)

site_order <- Fish_long_rel %>%
  distinct(Site) %>%
  mutate(
    Site_chr = as.character(Site),
    Prefix = substr(Site_chr, 1, 1),  # e.g., "O", "P"
    SampleType = stringr::str_extract(Site_chr, "_[A-Z]_?"),  
    Suffix = substr(Site_chr, nchar(Site_chr), nchar(Site_chr))  
  ) %>%
  arrange(
    factor(Prefix, levels = c("O", "P", "B", "S", "N")),
    SampleType,
    Suffix
  ) %>%
  pull(Site)

library(dplyr)
library(stringr)

site_order <- Fish_long_rel %>%
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



Fish_long_rel$Site <- factor(Fish_long_rel$Site, levels = site_order)
Fish_long_rel$Group <- factor(Fish_long_rel$Group, levels = c("O", "P", "B", "S"))

ggplot(Fish_long_rel, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(
    x = "Site",
    y = "Relative Abundance",
    fill = "Family"
  )

####SEPARATE FOR ACTIVE AND PASSIVE####
#to make it so active is combined change to "_L_|_S_" ~ "Active"
Fish_long_rel <- Fish_long_rel %>%
  mutate(SamplingType = case_when(
    str_detect(Site, "_M_") ~ "Passive",
    str_detect(Site, "_L_") ~ "Wilderlab",
    str_detect(Site, "_S_") ~ "Sterivex",
    TRUE ~ "Other"
  ))
table(Fish_long_rel$SamplingType)

####WILDERLAB####
wilder_data <- Fish_long_rel %>%
  filter(SamplingType == "Wilderlab") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Fish_long_rel %>%
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
    size = 9,
    vjust = 0
  ) +
geom_text(
    data = label_data,
    aes(x = Site, y = 14, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 9,
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
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Sample",
    y = "Relative Read Abundance",
    fill = "Family"
  )
p_wilder


##legend only
library(cowplot)

# Extract just the legend
legend <- get_legend(
  p_wilder + theme(legend.position = "right")
)

# Show only the legend
plot(legend)

####STERIVEX####
steri_data <- Fish_long_rel %>%
  filter(SamplingType == "Sterivex") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Fish_long_rel %>%
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
    size = 9,
    vjust = 0 ) +
  geom_text(
    data = label_data,
    aes(x = Site, y = 11, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 9,
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
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Sample",
    y = "Relative Read Abundance",
    fill = "Family"
  )
p_steri




####COMBINED ACTIVE####
Fish_long_rel <- Fish_long_rel %>%
  mutate(SamplingType = case_when(
    str_detect(Site, "_M_") ~ "Passive",
    str_detect(Site, "_L_|_S_") ~ "Active",
    TRUE ~ "Other"
  ))
table(Fish_long_rel$SamplingType)


active_data <- Fish_long_rel %>%
  filter(SamplingType == "Active") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Fish_long_rel %>%
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

total_abundance_per_site <- active_data %>%
  group_by(Site) %>%
  summarise(TotalAbundance = sum(RelAbundance))

# Plot
p_active <- ggplot(active_data, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
 geom_text(
    data = label_data %>% left_join(total_abundance_per_site, by = "Site"),
    aes(x = Site, y = 1.05 * TotalAbundance, label = FamilyCount),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
 geom_text(
    data = label_data,
    aes(x = Site, y = 12.5, label = SpeciesCount),
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
    y = "Relative Abundance",
    fill = "Family"
  )
p_active



passive_data <- Fish_long_rel %>%
  filter(SamplingType == "Passive") %>%
  mutate(
    Site = factor(Site),
    Group = factor(Group)
  )

family_counts <- Fish_long_rel %>%
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
total_abundance_per_site <- passive_data %>%
  group_by(Site) %>%
  summarise(TotalAbundance = sum(RelAbundance))

# Plot
p_passive <- ggplot(passive_data, aes(x = Site, y = RelAbundance, fill = ScientificName)) +
  geom_col(position = "stack") +
 geom_text(
    data = label_data %>% left_join(total_abundance_per_site, by = "Site"),
    aes(x = Site, y = 1.05 * TotalAbundance, label = FamilyCount),
    inherit.aes = FALSE,
    size = 9,
    vjust = 0
  ) +
 geom_text(
    data = label_data,
    aes(x = Site, y = 5, label = SpeciesCount),
    inherit.aes = FALSE,
    size = 9,
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
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Sample",
    y = "Relative Read Abundance",
    fill = "Family"
  )
p_passive

####USING PRIMERS####
Primers<-Fish

Primers[is.na(Primers)] <- 0
Primers[Primers == ""] <- 0

Primers[, -1] <- lapply(Primers[, -1], as.numeric)

print(Primers)

primers_12S <- Primers$Species[Primers$LM == 1]
primers_12S_ecoprimers <- Primers$Species[Primers$RV == 1]
primers_COI <- Primers$Species[Primers$CI == 1]
primers_16S <- Primers$Species[Primers$WV == 1]

venn_data <- list(
  '12S' = primers_12S,
  '12S ecoprimers' = primers_12S_ecoprimers,
  'COI' = primers_COI,
  '16S' = primers_16S
)

install.packages("ggVennDiagram")
library(ggVennDiagram)

####Gradient####
ggVennDiagram(venn_data, label = "count") +
  theme_void() +                         
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA) )
ggplot2::theme_minimal()

library(dplyr)

# Get all unique species detected by any primer
all_species <- unique(unlist(venn_data))

# Create binary presence/absence matrix
species_matrix <- data.frame(
  Species = all_species,
  RV = as.integer(all_species %in% venn_data$`12S ecoprimers`),
  WV = as.integer(all_species %in% venn_data$`16S`),
  CI = as.integer(all_species %in% venn_data$`COI`),
  LM = as.integer(all_species %in% venn_data$`12S`)
)

# Create a combination label 
species_matrix$Combination <- apply(species_matrix[, -1], 1, function(x) {
  present <- names(which(x == 1))
  if (length(present) == 0) {
    "None"
  } else {
    paste(present, collapse = "_")
  }
})

# Summary table: one row per unique combination
species_by_combination <- species_matrix %>%
  group_by(Combination) %>%
  summarise(
    Species = paste(Species, collapse = ", "),
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

# Print full table
print(species_by_combination)
####Jaccard####
presence_matrix <- Primers[, -1]  
presence_matrix[is.na(presence_matrix)] <- 0

t_matrix <- t(presence_matrix)
library(vegan)

# Calculate Jaccard distance (dissimilarity)
jaccard_dist <- vegdist(t_matrix, method = "jaccard", binary = TRUE)
jaccard_similarity <- 1 - as.matrix(jaccard_dist)

# Round and print
print(round(jaccard_similarity, 2))


library(ggplot2)
library(tidyr)
library(dplyr)

long_format <- species_matrix %>%
  pivot_longer(cols = RV:LM, names_to = "Primer", values_to = "Detected") %>%
  filter(Detected == 1)

species_order <- long_format %>%
  group_by(Species) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(Species)

long_format$Species <- factor(long_format$Species, levels = species_order)

# Plot
ggplot(long_format, aes(x = Primer, y = Species, fill = Primer)) +
  geom_tile(color = "white") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Species Detection by Primer",
    x = "Primer",
    y = "Species"
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

library(dplyr)
####USING GENERAL SHEET NMDS####
library(tidyverse)
library(ggVennDiagram)
install.packages("vegan")
library(vegan)
library(stringr)
####USING THE AGGREGATED ANALYSIS####
Fish_corrected <- Full_Fish_Dataset %>%
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

# Standardise each row (species) to 0–1
Fish_rel <- Fish_corrected %>%
  rowwise() %>%
  mutate(across(-ScientificName, ~ {
    vals <- c_across(-ScientificName)
    rng <- range(vals, na.rm = TRUE)
    if (diff(rng) == 0) 0 else (. - rng[1]) / diff(rng)
  })) %>%
  ungroup()

####RETRY####
rel_matrix <- Fish_rel %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame()
rel_matrix_t <- t(rel_matrix)

rel_matrix_t <- rel_matrix_t[rowSums(rel_matrix_t) > 0, ]
set.seed(123)
nmds <- metaMDS(rel_matrix_t, distance = "bray", k = 2, trymax = 100)
# Extract NMDS site scores
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Site <- rownames(site_scores)

# Clean names and extract site groups (e.g., O_L, P_S)
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

site_scores$NMDS1 <- -site_scores$NMDS1
ggplot(site_scores, aes(x = NMDS1, y = NMDS2,
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
metaprobe_cols <- names(Fish_rel)[str_detect(names(Fish_rel), "_M")]
df_meta <- Fish_rel %>%
  select(ScientificName, all_of(metaprobe_cols))

meta_matrix <- df_meta %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame() %>%
  t()  

meta_matrix <- meta_matrix[rowSums(meta_matrix) > 0, ]

set.seed(123)
nmds_meta <- metaMDS(meta_matrix, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores
meta_scores <- as.data.frame(scores(nmds_meta, display = "sites"))
meta_scores$Site <- rownames(meta_scores)

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
active_cols <- names(Fish_rel)[str_detect(names(Fish_rel), "_S|_L")]

df_active <- Fish_rel %>%
  select(ScientificName, all_of(active_cols))

active_matrix <- df_active %>%
  column_to_rownames("ScientificName") %>%
  as.data.frame() %>%
  t()  
active_matrix <- active_matrix[rowSums(active_matrix) > 0, ]

set.seed(123)
nmds_active <- metaMDS(active_matrix, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores
active_scores <- as.data.frame(scores(nmds_active, display = "sites"))
active_scores$Site <- rownames(active_scores)
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
# Sample names (from rel_matrix_t rownames)
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

# Identify which site(s) each species is associated with
# If 1 indicates presence and 0 absence
significant_results <- significant_results %>%
  mutate(Associated_Sites = paste0(
    ifelse(s.Bordjiesrif == 1, "Bordjiesrif, ", ""),
    ifelse(s.Olifants == 1, "Olifants, ", ""),
    ifelse(s.Platboom == 1, "Platboom, ", ""),
    ifelse(s.Smitswinkel == 1, "Smitswinkel, ", "")
  )) %>%
  mutate(Associated_Sites = gsub(", $", "", Associated_Sites)) # remove trailing comma

# Select and arrange columns
summary_table <- significant_results %>%
    select(Species, Associated_Sites, stat, p.value) %>%
    arrange(Associated_Sites, p.value)

# Make a nice table using gt
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

gtsave(my_gt_table, "indicator_species_table.png")



####REWORKED####
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(ggplot2)

# Filter to only indicator species
indicator_species <- summary_table$Species   

indicators_only <- rel_matrix_t[, colnames(rel_matrix_t) %in% indicator_species]

# Remove species with all zero counts
indicators_only <- indicators_only[, colSums(indicators_only) > 0]

# Build the OTU clustering tree
otu_dist <- vegdist(t(indicators_only), method = "bray")
otu_hclust <- hclust(otu_dist, method = "average")
otu_tree <- as.phylo(otu_hclust)

# Prepare metadata 
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
tip_metadata_long <- tip_metadata_long %>%
  left_join(data.frame(label = otu_tree$tip.label, y = 1:length(otu_tree$tip.label)), by = "label") %>%
  group_by(label) %>%
  mutate(circle_offset = seq_along(Associated_Sites) - 1)  
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

# Extract tip coordinates
tip_coords <- tree_base$data %>%
  filter(isTip) %>%
  select(label, x, y)

# Join metadata with coordinates
plot_data <- tip_metadata_long %>%
  left_join(tip_coords, by = "label") %>%
  mutate(x = x + 0.02 + 0.02 * stack)   

# Plot
tree_plot <- tree_base +
  # tip labels placed farthest out
  geom_tiplab(aes(label = label),
              size = 4,
              align = TRUE,
              offset = 0.07,   # move names outside
              linesize = 0.2) +
  # indicator circles close to tips
  geom_point(data = plot_data,
             aes(x = x, y = y, color = Associated_Sites),
             size = 4) +
  scale_color_manual(values = combined_colors, name = "Indicator Site(s)") +
  theme_void() +
  theme(legend.position = "right")

tree_plot

####lets see
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

# Indicator species 
indicator_species <- summary_table$Species

# All species matrix
all_species_matrix <- rel_matrix_t

# Remove species with all zero counts
all_species_matrix <- all_species_matrix[, colSums(all_species_matrix) > 0]

# Build OTU clustering tree 
otu_dist <- vegdist(t(all_species_matrix), method = "bray")
otu_hclust <- hclust(otu_dist, method = "average")
otu_tree <- as.phylo(otu_hclust)

# Metadata for all species
tip_metadata <- data.frame(
  label = otu_tree$tip.label,
  is_indicator = ifelse(otu_tree$tip.label %in% indicator_species, TRUE, FALSE),
  stringsAsFactors = FALSE
)

# For indicator species, attach their sites
tip_metadata <- tip_metadata %>%
  left_join(
    significant_results %>%
      select(Species, Associated_Sites) %>%
      rename(label = Species),
    by = "label"
  )

tip_metadata$Associated_Sites[is.na(tip_metadata$Associated_Sites)] <- "None"

# Split multi-site species into separate rows for stacking
tip_metadata_long <- tip_metadata %>%
  separate_rows(Associated_Sites, sep = ",\\s*") %>%
  group_by(label) %>%
  mutate(stack = row_number() - 1) %>%
  ungroup()

# ---- Custom colors ----
manual_colors <- c(
  "Bordjiesrif" = "#CC79A7",
  "Olifants" = "blue",
  "Platboom" = "orange",
  "Smitswinkel" = "green",
  "Olifants, Platboom"= "yellow",
  "Bordjiesrif, Smitswinkel"="red"
)

all_categories <- unique(tip_metadata_long$Associated_Sites)
auto_cats <- setdiff(all_categories, names(manual_colors))
auto_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(auto_cats)),
  auto_cats
)
combined_colors <- c(manual_colors, auto_colors)
combined_colors["None"] <- "black"

tree_base <- ggtree(otu_tree, layout = "circular")

tip_coords <- tree_base$data %>%
  filter(isTip) %>%
  select(label, x, y)

# Join metadata with coordinates 
plot_data <- tip_metadata_long %>%
  left_join(tip_coords, by = "label") %>%
  mutate(x = x + 0.02 + 0.02 * stack)

# Plot tree 

tree_plot <- tree_base +
  geom_tiplab(aes(label = label),
              size = 6,
              align = TRUE,
              offset = 0.07,
              linesize = 0.2) +
  geom_point(data = plot_data,
             aes(x = x, y = y, color = Associated_Sites),
             size = 6) +
  scale_color_manual(values = combined_colors, name = "Indicator Site(s)") +
  theme_void() +
  theme(legend.position ="none")

tree_plot

ggsave(
  filename = "tree_plot.png",
  plot = tree_plot,
  width = 12,           
  height = 12,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)
