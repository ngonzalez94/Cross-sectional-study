#Load packages 
library(readxl)
library(tidyverse)
library(vegan)  # for vegdist and other ecological analyses
install.packages("ape")
library(ape)    # for pcoa

# Define the file path
file_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken/Relative abundance PCoA.xlsx"

# Load the bracken data and convert to tibble
dbf <- read_excel(file_path, sheet=2) %>% as_tibble()
dkf <- read_excel(file_path, sheet=4) %>% as_tibble()
drm <- read_excel(file_path, sheet=8) %>% as_tibble()
dwm <- read_excel(file_path, sheet=6) %>% as_tibble()

# Function to prepare data for PCoA
prepare_pcoa_data <- function(df, sample_prefix) {
  df <- as.data.frame(df)
  
  # Get all sample columns (both T1 and T2)
  sample_cols <- grep(paste0(sample_prefix, ".*_T[12]$"), colnames(df), value = TRUE)
  
  # Extract abundance data and transpose (samples as rows, taxa as columns)
  abundance_matrix <- df %>%
    select(all_of(sample_cols)) %>%
    t() %>%
    as.data.frame()
  
  # Set column names as taxa names
  colnames(abundance_matrix) <- df[[1]]
  
  # Add sample information
  abundance_matrix$Sample <- rownames(abundance_matrix)
  abundance_matrix$Time <- ifelse(grepl("_T1$", rownames(abundance_matrix)), "Time1", "Time2")
  abundance_matrix$Group <- sample_prefix
  
  return(abundance_matrix)
}

# Prepare data for all samples
bf_data <- prepare_pcoa_data(dbf, "B")
kf_data <- prepare_pcoa_data(dkf, "K")
rm_data <- prepare_pcoa_data(drm, "RM")
wm_data <- prepare_pcoa_data(dwm, "WM")

# Combine all data
all_data <- bind_rows(bf_data, kf_data, rm_data, wm_data)

# Separate metadata from abundance data
metadata <- all_data %>%
  select(Sample, Time, Group)

abundance_data <- all_data %>%
  select(-Sample, -Time, -Group)

# Replace NAs with 0
abundance_data[is.na(abundance_data)] <- 0

# Calculate Bray-Curtis dissimilarity
bray_dist <- vegdist(abundance_data, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_dist)

# Extract coordinates
pcoa_coords <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")

# Add metadata
pcoa_coords <- cbind(pcoa_coords, metadata)

# Calculate variance explained
var_explained <- pcoa_result$values$Relative_eig
pc1_var <- round(var_explained[1] * 100, 2)
pc2_var <- round(var_explained[2] * 100, 2)

# Create PCoA plot
ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, color = Group, shape = Time)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95, 
               linetype = 2, linewidth = 0.8)+
  scale_color_manual(values = c("B" = "#E41A1C", "K" = "#377EB8", 
                                "RM" = "#4DAF4A", "WM" = "#984EA3"),
                     labels = c("B" = "CMR-fed calves' faeces", "K" = "WM-fed calves' faeces", 
                                "RM" = "CMR", "WM" = "WM")) +
  scale_shape_manual(values = c("Time1" = 16, "Time2" = 17),
                     labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  theme_minimal() +
  labs(
    x = paste0("PCoA1 (", pc1_var, "%)"),
    y = paste0("PCoA2 (", pc2_var, "%)"),
    color = "Sample Type",
    shape = "Time Point"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size =13),
    legend.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )


# After running your PCoA analysis...

# 1. Get the variance explained by each axis
var_explained <- pcoa_result$values$Relative_eig
cat("Variance explained by axes:\n")
cat("PCoA1:", round(var_explained[1] * 100, 2), "%\n")
cat("PCoA2:", round(var_explained[2] * 100, 2), "%\n")
cat("Total (PC1+PC2):", round(sum(var_explained[1:2]) * 100, 2), "%\n")

# 2. Correlate taxa with PCoA axes to see what drives separation
# Get PCoA scores
pcoa_scores <- pcoa_result$vectors[, 1:2]

# Calculate correlations between taxa abundances and PCoA axes
taxa_correlations <- cor(abundance_data, pcoa_scores)
colnames(taxa_correlations) <- c("PCoA1", "PCoA2")

# Get top taxa correlated with each axis
top_n <- 10

cat("\nTop taxa correlated with PCoA1:\n")
taxa_pcoa1 <- taxa_correlations[order(abs(taxa_correlations[, "PCoA1"]), decreasing = TRUE), ]
print(head(taxa_pcoa1, top_n))

cat("\nTop taxa correlated with PCoA2:\n")
taxa_pcoa2 <- taxa_correlations[order(abs(taxa_correlations[, "PCoA2"]), decreasing = TRUE), ]
print(head(taxa_pcoa2, top_n))

# 3. Statistical test: PERMANOVA to test if groups differ significantly
library(vegan)

# Test effect of Group
permanova_group <- adonis2(bray_dist ~ Group, data = metadata, permutations = 999)
print("PERMANOVA - Effect of Sample Type:")
print(permanova_group)

# Test effect of Time
permanova_time <- adonis2(bray_dist ~ Time, data = metadata, permutations = 999)
print("PERMANOVA - Effect of Time:")
print(permanova_time)

# Test both factors
permanova_both <- adonis2(bray_dist ~ Group * Time, data = metadata, permutations = 999)
print("PERMANOVA - Effect of Sample Type and Time:")
print(permanova_both)

# 4. Test for homogeneity of dispersion (important for PERMANOVA interpretation)
# Calculate dispersions
disp_group <- betadisper(bray_dist, metadata$Group)
disp_time <- betadisper(bray_dist, metadata$Time)

# Test for significant differences in dispersion
print("Homogeneity of dispersion - Group:")
print(permutest(disp_group))

print("Homogeneity of dispersion - Time:")
print(permutest(disp_time))



# ============================================================
# NMDS Analysis (add this after your PCoA code)
# ============================================================

# Perform NMDS
set.seed(123)  # for reproducibility
nmds_result <- metaMDS(abundance_data, distance = "bray", k = 2, trymax = 100)

# Extract stress value
stress_value <- round(nmds_result$stress, 3)
cat("NMDS Stress:", stress_value, "\n")
cat("Stress interpretation: <0.05 excellent, <0.1 good, <0.2 acceptable, >0.2 poor\n")

# Extract NMDS coordinates
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
colnames(nmds_coords) <- c("NMDS1", "NMDS2")

# Add metadata
nmds_coords <- cbind(nmds_coords, metadata)

# Create NMDS plot
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Group)) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95, 
               linetype = 2, linewidth = 0.8) +
  geom_point(aes(shape = Time), size = 4, alpha = 0.8) +
  scale_color_manual(values = c("B" = "#E41A1C", "K" = "#377EB8", 
                                "RM" = "#4DAF4A", "WM" = "#984EA3"),
                     labels = c("B" = "BF", "K" = "KF", 
                                "RM" = "RM", "WM" = "WM")) +
  scale_shape_manual(values = c("Time1" = 16, "Time2" = 17),
                     labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  theme_minimal() +
  labs(
    title = paste0("NMDS of Microbial Communities (Stress: ", stress_value, ")"),
    x = "NMDS1",
    y = "NMDS2",
    color = "Sample Type",
    shape = "Time Point"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )
