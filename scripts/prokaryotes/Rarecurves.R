library(tidyverse)
library(vegan)

####Rarecurves by site####
# Convert to matrix
otu_matrix <- as.matrix(Unedited_Bacteria[, -1])
rownames(otu_matrix) <- Unedited_Bacteria$ScientificName

# Transpose for rarecurve: samples in rows, taxa (OTUs) in columns
otu_t <- t(otu_matrix)

# Remove lab blanks (assuming names start with "N")
otu_t <- otu_t[!grepl("^N", rownames(otu_t)), ]
# Extract first character as site code
site_codes <- substr(rownames(otu_t), 1, 1)

# Map short codes to full site names
site_map <- c(
  O = "Olifantsbos",
  P = "Platboom",
  B = "Bordjiesrif",
  S = "Smitswinkel"
)

# Get full site names
site_names <- site_map[site_codes]

# Create site-wise pooled abundance data (sum OTUs per site)
pooled_by_site <- rowsum(otu_t, group = site_names)

# Assign site colors
site_colors <- setNames(rainbow(nrow(pooled_by_site)), rownames(pooled_by_site))

# Plot rarefaction curves per site
# Suppress direct labels on lines
dev.off() 

rarecurve(
  pooled_by_site,
  step = 1000,
  sample = min(rowSums(pooled_by_site)),
  col = site_colors,
  label = FALSE,  # disable messy labels
  xlab = "Reads (pooled per site)",
  ylab = "OTU richness",
  main = "Rarefaction Curves (Pooled by Site)",
  lwd = 2
)

# Add a clear, non-overlapping legend
legend(
  "bottomright",                        # or "topright" if clearer
  legend = names(site_colors),         # Site names
  col = site_colors,                   # Same color mapping
  lty = 1,                             # Line type
  lwd = 2,                             # Line width to match plot
  bty = "n",
  y.intersp = 0.2 
)

