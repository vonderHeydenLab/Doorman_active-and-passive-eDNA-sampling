####sort data####
library(dplyr)
library(tidyverse)
library(ggVennDiagram)
library(VennDiagram)

####Olifants####
Olifants <- Fish %>%
  mutate(
    O_M = rowSums(select(., starts_with("O_M"))) > 0,
    O_S = rowSums(select(., starts_with("O_S"))) > 0,
    O_L = rowSums(select(., starts_with("O_L"))) > 0
  ) %>%
  select(ScientificName, O_M, O_S, O_L)

# Create lists of species per site
Olifants_list <- list(
  Metaprobe = Olifants$ScientificName[Olifants$O_M],
  Sterivex = Olifants$ScientificName[Olifants$O_S],
  Wilderlab = Olifants$ScientificName[Olifants$O_L]
)

ggVennDiagram(Olifants_list, label = "count") +
  scale_fill_gradient(low = "white", high = "dodgerblue") 

####Platboom####
library(VennDiagram)
Platboom <- Fish %>%
  mutate(
    P_M = rowSums(select(., starts_with("P_M"))) > 0,
    P_S = rowSums(select(., starts_with("P_S"))) > 0
  ) %>%
  select(ScientificName, P_M, P_S)

Platboom_list <- list(
  Metaprobe = Platboom$ScientificName[Platboom$P_M],
  Sterivex  = Platboom$ScientificName[Platboom$P_S]
)

venn.plot <- venn.diagram(
  x = Platboom_list,
  category.names = c("Metaprobe", "Sterivex"),
  filename = NULL,  
  output = TRUE,
  
  # Aesthetic settings
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.5,            
  cat.cex = 1.2,        
  cat.pos = c(180, 0),  
  cat.dist = c(0.05,0.05),
  margin = 0.2
)

grid.newpage()
grid.draw(venn.plot)

####Bordjiesrif####
Bordjiesrif <- Fish %>%
  mutate(
    B_M = rowSums(select(., starts_with("B_M"))) > 0,
    B_S = rowSums(select(., starts_with("B_S"))) > 0
  ) %>%
  select(ScientificName, B_M, B_S)

Bordjiesrif_list <- list(
  Metaprobe = Bordjiesrif$ScientificName[Bordjiesrif$B_M],
  Sterivex  = Bordjiesrif$ScientificName[Bordjiesrif$B_S]
)

venn.plot <- venn.diagram(
  x = Bordjiesrif_list,
  category.names = c("Metaprobe", "Sterivex"),
  filename = NULL,  
  output = TRUE,
  
  # Aesthetic settings
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.5,            
  cat.cex = 1.2,        
  cat.pos = c(180, 0),  
  cat.dist = c(0.08,0.08),
  margin = 0.2
)

grid.newpage()
grid.draw(venn.plot)

####Smitswinkel####
Smits <- Fish %>%
  mutate(
    S_M = rowSums(select(., starts_with("S_M"))) > 0,
    S_S = rowSums(select(., starts_with("S_S"))) > 0,
    S_L = rowSums(select(., starts_with("S_L"))) > 0
  ) %>%
  select(ScientificName, S_M, S_S, S_L)

# Create lists of species per site
Smits_list <- list(
  Metaprobe = Smits$ScientificName[Smits$S_M],
  Sterivex = Smits$ScientificName[Smits$S_S],
  Wilderlab = Smits$ScientificName[Smits$S_L]
)

ggVennDiagram(Smits_list, label = "count") +
  scale_fill_gradient(low = "white", high = "dodgerblue") 


####Combined image####
library(VennDiagram)
library(gridExtra)
library(grid)

# Create all four venn.diagram grobs

Olifants_venn <- venn.diagram(
  x = Olifants_list,
  category.names = c("Metaprobe", "Sterivex", "Wilderlab"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue", "lightgreen"),
  alpha = 0.5,
  cex = 1.2, 
  cat.pos = c(0, 120, 240),   # Adjust angles for each label
  cat.dist = c(0.08, 0.08, 0.08), 
  cat.cex = 1, 
  margin = 0.1
)

Platboom_venn <- venn.diagram(
  x = Platboom_list,
  category.names = c("Metaprobe", "Sterivex"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1, 
  cat.dist = c(0.08,0.08),
  margin = 0.1
)

Bordjiesrif_venn <- venn.diagram(
  x = Bordjiesrif_list,
  category.names = c("Metaprobe", "Sterivex"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.2, 
  cat.dist = c(0.08,0.08),
  cat.cex = 1,
  margin = 0.1
)

Smits_venn <- venn.diagram(
  x = Smits_list,
  category.names = c("Metaprobe", "Wilderlab", "Sterivex"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue", "lightgreen"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1,
  cat.pos = c(20, 120, 240),   # Adjust angles for each label
  cat.dist = c(0.05, 0.2, 0.1), 
  margin = 0.1
)
grid.newpage()
grid.draw(Smits_venn)

# Arrange all 4 plots in a 2x2 grid
grid.newpage()
grid.arrange(
  grobs = list(Olifants_venn, Platboom_venn, Bordjiesrif_venn, Smits_venn),
  ncol = 2,
  nrow = 2
)

library(grid)
library(gridExtra)
library(VennDiagram)

# Your Venn grobs (e.g., Olifants_venn, Platboom_venn, etc.) created previously

# Create text grobs for labels
label_1 <- textGrob("Olifants", gp = gpar(fontsize = 14),vjust = 0.4)
label_2 <- textGrob("Platboom", gp = gpar(fontsize = 14),vjust=0.4)
label_3 <- textGrob("Bordjiesrif", gp = gpar(fontsize = 14), vjust=0.3)
label_4 <- textGrob("Smitswinkel", gp = gpar(fontsize = 14),vjust=0.3)

# Combine each venn with its label using arrangeGrob vertically
plot1 <- arrangeGrob(Olifants_venn, label_1, ncol = 1, heights = c(4, 0.15))
plot2 <- arrangeGrob(Platboom_venn, label_2, ncol = 1, heights = c(4, 0.15))
plot3 <- arrangeGrob(Bordjiesrif_venn, label_3, ncol = 1, heights = c(4, 0.25))
plot4 <- arrangeGrob(Smits_venn, label_4, ncol = 1, heights = c(4, 0.25))

# Arrange all four labeled plots in 2x2 grid
grid.newpage()
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)


####Comparing Active Samples####
Active_Oli <- list(
  Sterivex = Olifants$ScientificName[Olifants$O_S],
  Wilderlab = Olifants$ScientificName[Olifants$O_L]
)

OliAct <- venn.diagram(
  x = Active_Oli,
  category.names = c("Sterivex", "Wilderlab"),
  filename = NULL,  
  output = TRUE,
  
  # Aesthetic settings
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.5,            
  cat.cex = 1.2,        
  cat.pos = c(220,-220),  
  cat.dist = c(0.05,0.05),
  margin = 0.2
)

grid.newpage()
grid.draw(OliAct)


Active_Smits <- list(
  Sterivex = Smits$ScientificName[Smits$S_S],
  Wilderlab = Smits$ScientificName[Smits$S_L]
)

SmitsAct <- venn.diagram(
  x = Active_Smits,
  category.names = c("Sterivex", "Wilderlab"),
  filename = NULL,  
  output = TRUE,
  
  # Aesthetic settings
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.5,            
  cat.cex = 1.2,        
  cat.pos = c(220,-220),  
  cat.dist = c(0.05,0.05),
  margin = 0.2
)

grid.newpage()
grid.draw(SmitsAct)

####Active Time####
Day1 <- Fish %>%
  mutate(
    SteriOli1 = rowSums(select(., "O_S_1_1","O_S_2_1","O_S_3_1")) > 0,
    Plat1 = rowSums(select(.,"P_S_1_1","P_S_2_1","P_S_3_1" )) > 0,
    Borj1 = rowSums(select(.,"B_S_1_1","B_S_2_1","B_S_3_1" )) > 0,
    SteriSmit1 = rowSums(select(.,"S_S_1_1","S_S_2_1","S_S_3_1" )) > 0,
    WilderSmit1 = rowSums(select(.,"S_L_1_1","S_L_2_1","S_L_3_1" )) > 0,
    WilderOli1 = rowSums(select(.,"O_L_1_1", "O_L_2_1","O_L_3_1")) > 0
  ) %>%
  select(ScientificName, SteriOli1, Plat1, Borj1, SteriSmit1, WilderSmit1, WilderOli1)

Day2<- Fish %>%
  mutate(
    SteriOli2 = rowSums(select(., "O_S_1_2","O_S_2_2","O_S_3_2")) > 0,
    Plat2 = rowSums(select(.,"P_S_1_2","P_S_2_2","P_S_3_2" )) > 0,
    Borj2 = rowSums(select(.,"B_S_1_2","B_S_2_2","B_S_3_2" )) > 0,
    SteriSmit2 = rowSums(select(.,"S_S_1_2","S_S_2_2","S_S_3_2" )) > 0,
    WilderSmit2 = rowSums(select(.,"S_L_1_2","S_L_2_2","S_L_3_2" )) > 0,
    WilderOli2 = rowSums(select(.,"O_L_1_2", "O_L_2_2","O_L_3_2")) > 0
  ) %>%
  select(ScientificName, SteriOli2, Plat2, Borj2, SteriSmit2, WilderSmit2, WilderOli2)

####Comparing Day 1 and Day 2 ####
Day1$Detected <- Day1$SteriOli1 | Day1$Plat1 | Day1$Borj1 | 
  Day1$SteriSmit1 | Day1$WilderSmit1 | Day1$WilderOli1

Day2$Detected <- Day2$SteriOli2 | Day2$Plat2 | Day2$Borj2 | 
  Day2$SteriSmit2 | Day2$WilderSmit2 | Day2$WilderOli2

venn_days <- list(
  Day1 = Day1$ScientificName[Day1$Detected],
  Day2 = Day2$ScientificName[Day2$Detected]
)

ggVennDiagram(venn_days, label = "count") +
  theme_void() +
  ggtitle("Species Detected on Day 1 vs Day 2")


species_day1 <- Day1$ScientificName[Day1$Detected]
species_day2 <- Day2$ScientificName[Day2$Detected]

unique_to_day1 <- setdiff(species_day1, species_day2)
unique_to_day2 <- setdiff(species_day2, species_day1)
shared_species <- intersect(species_day1, species_day2)

length(unique_to_day1)      # How many unique to Day 1
length(unique_to_day2)      # How many unique to Day 2
length(shared_species)      # How many shared

# Print some example species
print(unique_to_day1)
print(unique_to_day2)
print(shared_species)


combined_species <- data.frame(
  ScientificName = c(unique_to_day1, unique_to_day2, shared_species),
  Category = c(
    rep("Day1 only", length(unique_to_day1)),
    rep("Day2 only", length(unique_to_day2)),
    rep("Shared", length(shared_species))
  )
)
write.csv(combined_species, "combined_species_day1_day2.csv", row.names = FALSE)
