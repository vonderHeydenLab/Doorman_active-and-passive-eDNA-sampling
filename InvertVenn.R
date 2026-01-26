library(dplyr)
library(tidyverse)
library(ggVennDiagram)
library(VennDiagram)

Olifants <- Invert %>%
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

Platboom <- Invert %>%
  mutate(
    P_M = rowSums(select(., starts_with("P_M"))) > 0,
    P_S = rowSums(select(., starts_with("P_S"))) > 0
  ) %>%
  select(ScientificName, P_M, P_S)

Platboom_list <- list(
  Metaprobe = Platboom$ScientificName[Platboom$P_M],
  Sterivex  = Platboom$ScientificName[Platboom$P_S]
)

Bordjiesrif <- Invert %>%
  mutate(
    B_M = rowSums(select(., starts_with("B_M"))) > 0,
    B_S = rowSums(select(., starts_with("B_S"))) > 0
  ) %>%
  select(ScientificName, B_M, B_S)

Bordjiesrif_list <- list(
  Metaprobe = Bordjiesrif$ScientificName[Bordjiesrif$B_M],
  Sterivex  = Bordjiesrif$ScientificName[Bordjiesrif$B_S]
)

Smits <- Invert %>%
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
  cat.pos = c(0, 95, 240),   # Adjust angles for each label
  cat.dist = c(0.08, 0.08, 0.08), 
  cat.cex = 1, 
  margin = 0.1
)
grid.newpage()
grid.draw(Olifants_venn)

Platboom_venn <- venn.diagram(
  x = Platboom_list,
  category.names = c("Metaprobe", "Sterivex"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1, 
  cat.pos = c(-90, 90),
  cat.dist = c(0.08,0.08),
  margin = 0.1
)
grid.newpage()
grid.draw(Platboom_venn)


Bordjiesrif_venn <- venn.diagram(
  x = Bordjiesrif_list,
  category.names = c("Metaprobe", "Sterivex"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue"),
  alpha = 0.5,
  cex = 1.2, 
  cat.pos = c(-90, 90),
  cat.dist = c(0.08,0.08),
  cat.cex = 1,
  margin = 0.1
)
grid.newpage()
grid.draw(Bordjiesrif_venn)

Smits_venn <- venn.diagram(
  x = Smits_list,
  category.names = c("Metaprobe", "Wilderlab", "Sterivex"),
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "skyblue", "lightgreen"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1,
  cat.pos = c(0, 95, 240),   # Adjust angles for each label
  cat.dist = c(0.08, 0.08, 0.08), 
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
label_3 <- textGrob("Bordjiesrif", gp = gpar(fontsize = 14), vjust=0.05)
label_4 <- textGrob("Smitswinkel", gp = gpar(fontsize = 14),vjust=0.05)

# Combine each venn with its label using arrangeGrob vertically
plot1 <- arrangeGrob(Olifants_venn, label_1, ncol = 1, heights = c(4, 0.15))
plot2 <- arrangeGrob(Platboom_venn, label_2, ncol = 1, heights = c(4, 0.15))
plot3 <- arrangeGrob(Bordjiesrif_venn, label_3, ncol = 1, heights = c(4, 0.1))
plot4 <- arrangeGrob(Smits_venn, label_4, ncol = 1, heights = c(4, 0.1))

# Arrange all four labeled plots in 2x2 grid
grid.newpage()
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)



