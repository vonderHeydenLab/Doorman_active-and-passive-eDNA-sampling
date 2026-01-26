library(tidyverse)
library(vegan)

# Prepare OTU matrix 
otu_matrix <- as.matrix(Full_Fish_Dataset[, -1])
rownames(otu_matrix) <- Full_Fish_Dataset$ScientificName
otu_t <- t(otu_matrix)

# Remove lab blanks and all-zero samples 
otu_t <- otu_t[!grepl("^N", rownames(otu_t)), ]
otu_t <- otu_t[rowSums(otu_t) > 0, ]  

# Get site names 
site_codes <- substr(rownames(otu_t), 1, 1)
site_map <- c(O = "Olifantsbos", P = "Platboom", B = "Bordjiesrif", S = "Smitswinkel")
site_names <- site_map[site_codes]

# Rarefy each sample 
min_depth <- min(rowSums(otu_t))
cat("Minimum sequencing depth used for rarefaction:", min_depth, "\n")

otu_rarefied <- rrarefy(otu_t, sample = min_depth)

# Pool by site 
pooled_rarefied <- rowsum(otu_rarefied, group = site_names)

# Assign colors 
site_colors <- setNames(rainbow(length(unique(site_names))), unique(site_names))
sample_colors <- site_colors[site_names]

# Plot rarecurves for individual samples
rarecurve(otu_rarefied,
          step = 1000,
          col = sample_colors,
          label = FALSE,
          xlab = "Number of Sequences",
          ylab = "Number of OTUs",
          main = "Rarefaction Curves: Individual Rarefied Samples")

legend("bottomright", legend = names(site_colors), col = site_colors, lwd = 2, bty = "n")

#Plot rarecurves for pooled samples
rarecurve(pooled_rarefied,
          step = 1000,
          col = site_colors,
          label = FALSE,
          xlab = "Number of Sequences",
          ylab = "Number of OTUs",
          main = "Rarefaction Curves: Pooled Rarefied Samples by Site")
legend("bottomright", legend = names(site_colors), col = site_colors, lwd = 2, bty = "n",
       y.intersp = 0.2 )

