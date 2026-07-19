# Load required libraries
library(tidyverse)
library(viridis)

# Set working directory
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Long_reads/")

# Specify the files you want to include
files <- c("MD0117_ARGs_with_families.csv",
           "E06_ARGs_with_families.csv")

# Check which files exist
files_exist <- file.exists(files)
files <- files[files_exist]

print(paste("Found", length(files), "files:"))
print(files)

# Read and combine files
combined_data <- data.frame()

for (file in files) {
  # Read the CSV
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract sample ID from filename (e.g., "MD0117" from "MD0117_ARGs_with_families.csv")
  sample_id <- gsub("_ARGs_with_families\\.csv$", "", file)
  
  # Add sample ID column
  df$Sample <- sample_id
  
  # Combine
  combined_data <- rbind(combined_data, df)
}

print(paste("Total rows loaded:", nrow(combined_data)))

# Clean column names (remove extra spaces)
colnames(combined_data) <- trimws(colnames(combined_data))

# Print column names to verify
print("Column names:")
print(colnames(combined_data))

# Check if key columns exist
if (!"GENE" %in% colnames(combined_data)) {
  stop("'GENE' column not found. Available columns: ", paste(colnames(combined_data), collapse = ", "))
}
if (!"Family" %in% colnames(combined_data)) {
  stop("'Family' column not found. Available columns: ", paste(colnames(combined_data), collapse = ", "))
}

# Remove rows with missing values in key columns
combined_data <- combined_data %>%
  filter(!is.na(GENE), !is.na(Family), !is.na(Sample))

print(paste("Rows after removing NAs:", nrow(combined_data)))

# Remove viral families (ending in "viridae")
#combined_data <- combined_data %>%
#  filter(!str_detect(Family, "viridae$"))

#print(paste("Rows after removing viridae:", nrow(combined_data)))

# Check unique samples
print("Unique samples:")
print(unique(combined_data$Sample))

# Add custom labels for samples
combined_data <- combined_data %>%
  mutate(Sample_Label = case_when(
     Sample == "MD0117" ~ "Replacement calves faeces Farm2 Time 1 (Long reads)",
    Sample == "E06" ~ "Replacement calves faeces Farm1 Time 1 (Short reads)",
    TRUE ~ Sample  # Keep original name if not matched
  ))%>%
  mutate(Sample_Label = factor(Sample_Label, 
                               levels = c("Replacement calves faeces Farm2 Time 1 (Long reads)",
                                          "Replacement calves faeces Farm1 Time 1 (Short reads)")))

print("Sample labels:")
print(unique(combined_data$Sample_Label))

# Check unique families
print(paste("Number of unique families:", length(unique(combined_data$Family))))

# Check unique genes
print(paste("Number of unique genes:", length(unique(combined_data$GENE))))

# Count occurrences of each Gene-Family-Sample combination
heatmap_data <- combined_data %>%
  group_by(Sample_Label, Family, GENE) %>%
  summarise(Count = n(), .groups = "drop")

print(paste("Rows in heatmap_data:", nrow(heatmap_data)))
print(paste("Number of unique genes:", length(unique(heatmap_data$GENE))))
print(paste("Number of unique families:", length(unique(heatmap_data$Family))))

# Get top families by total count (to reduce clutter)
top_families <- heatmap_data %>%
  group_by(Family) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  head(20) %>%  # Keep top 20 families
  pull(Family)

# Filter to top families
heatmap_plot_data <- heatmap_data %>%
  filter(Family %in% top_families)

print(paste("Final plot data rows:", nrow(heatmap_plot_data)))
print(paste("Genes to plot:", length(unique(heatmap_plot_data$GENE))))

# Create the heatmap with individual gene names
p <- ggplot(heatmap_plot_data, aes(x = Family, y = GENE, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ Sample_Label, ncol = 3, scales = "free_x") +  # Use Sample_Label for custom titles
  scale_fill_viridis(option = "plasma", name = "Count", direction = -1) +
  theme_minimal() +
  labs(
    #title = "Antimicrobial Resistance Genes by Taxonomic Family",
   #subtitle = "Only showing families with ARGs in each sample",
    x = "Taxonomic Family",
    y = "Resistance Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

# Save the plot - adjust height based on number of genes
plot_height <- max(10, length(unique(heatmap_plot_data$GENE)) * 0.2)
ggsave("ARG_by_Taxa_Heatmap_F2_T1.png", plot = p, width = 16, height = plot_height, dpi = 300, limitsize = FALSE)
print(paste("Plot saved as ARG_by_Taxa_Heatmap.png with height:", plot_height))

print("Done!")
print(p)
