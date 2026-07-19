library(tidyverse)
library(readxl)
library(readr)
library(ggplot2)
library(dplyr)

# Load ggplot2
# Load your ARG read count data
arg_counts <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/16S_ARG/ARG_counts_sample_gene_noSNP.xlsx", sheet = 1)

# Extract categories (I'm deleting Thiopeptides and Tetracenomycin and Spiropurimidinetriones manually because they are NA)
arg_counts <- arg_counts %>%
  filter(!grepl("Biocides", Gene_Name) & (!grepl("Metals", Gene_Name) & (!grepl("16S", Gene_Name) & (!grepl("metal", Gene_Name)) & (!grepl("G07", Gene_Name) & (!grepl("Thiopeptides", Gene_Name) & (!grepl("Tetracenomycin", Gene_Name) & (!grepl("Spiropyrimidinetriones", Gene_Name))))))))

# Extract gene name from the 'Gene' column (5th part of the pipe-delimited name) So it deletes the rest
arg_counts <- arg_counts %>%
  mutate(
    Gene_Name = sapply(strsplit(Gene_Name, "\\|"), function(x) paste(x[c(1, 3)], collapse = "|"))
  )


# Load ARG gene lengths (e.g., from MEGARes FASTA processed via bioawk) This renames the columns into Gene name (everything from MEG to the actual name of the gene) and then the arg length)
arg_lengths <- read_tsv("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/16S_ARG/arg_lengths.tsv", col_names = c("Gene_Name", "ARG_Length"))

# Extract gene name from the 'Gene' column (5th part of the pipe-delimited name) So it deletes the rest
# Extract gene name from the 'Gene' column (5th part of the pipe-delimited name) So it deletes the rest (Inside the Gene name column we stay with the divisions one and three)
arg_lengths <- arg_lengths %>%
  mutate(
    Gene_Name = sapply(strsplit(Gene_Name, "\\|"), function(x) paste(x[c(1, 3)], collapse = "|"))
  )

# Load 16S counts
s16_data <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/16S_ARG/16S_counts.xlsx")

# Merge ARG read data with ARG lengths
arg_data <- left_join(arg_counts, arg_lengths, by = "Gene_Name", relationship = "many-to-many")


# Merge with 16S count data
merged <- left_join(arg_data, s16_data, by = "Sample")

# Constants
L_reads <- 150  # average read length (adjust if needed)
L_16S <- 1432   # 16S reference length as used per Greengenes

#arg_summary <- merged %>%
#  group_by(Sample, Gene_Name) %>%
#  summarise(
#    Total_Reads = sum(Read_Counts / ARG_Length, na.rm = TRUE),
#    Total_16S = unique(`16S_Count` / L_16S, na.rm = TRUE),
#    .groups = "drop"
#  ) %>%
#  mutate(ARG_per_16S = Total_Reads / Total_16S)

arg_summary <- merged %>%
  mutate(
    Class = sapply(strsplit(Gene_Name, "\\|"), function(x) x[2])
  ) %>%
  group_by(Sample, Class) %>%   # AGGREGATE BY CLASS HERE
  summarise(
    Total_Reads = sum(Read_Counts*L_reads/ ARG_Length, na.rm = TRUE),
    Total_16S = unique(`16S_Count`*L_reads / L_16S, na.rm = TRUE),
    ARG_per_16S = Total_Reads / Total_16S,
    .groups = "drop"
  )

#Assign categories manually
arg_summary <- arg_summary %>%
  mutate(Category = case_when(
    Sample %in% c("E06", "E02", "E10", "F01", "F04", "F12", "G04", "G09", "G12", "H04", "H08") ~ "Non-replacement calves' faeces",
    Sample %in% c("E03", "E07", "F05", "F09", "G01", "G05", "H01", "H05", "H09") ~ "Replacement milk",
    Sample %in% c("E01", "E05", "E09", "E12", "F03", "F11", "G03", "G08", "G11", "H03", "H07") ~ "Replacement calves' faeces",
    Sample %in% c("E04", "E08", "E11", "F02", "F06", "F10", "G02", "G06", "G10", "H06", "H10") ~ "Waste milk",
    TRUE ~ "Unknown"
  ))

# Example data frame: assume you have columns Sample, gene_name, Ratio, etc.

# Step 1: Extract Class from gene_name (assuming class is somewhere in gene_name string)
# For example, if gene_name is "MEGXX|aminoglycoside|something", extract "aminoglycoside"

#arg_summary <- arg_summary %>%
 # mutate(Class = sapply(strsplit(Class, "\\|"), function(x) x[2]))  # Adjust the index based on your data

arg_summary <- arg_summary %>%
  filter(Sample != "G07")


#log_scale
arg_summary <- arg_summary %>%
  mutate(log10_ARG_per_16S = log10(ARG_per_16S + 1))

##plot_logged
ggplot(arg_summary, aes(x = Category, y = Class, fill = Category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  scale_y_log10() +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  labs(
    y = "ARG per 16S rRNA gene (log10 scale)",
    x = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

#Permanently remove a group from the x-axis
arg_summary_filtered <- arg_summary %>%
  filter(Class != "Drug_and_biocide_resistance")


# Step 2: Now group by Class and create boxplot for Ratio (or whatever metric you have)
ggplot(arg_summary_filtered, aes(x = ARG_per_16S, y = Class)) +
  geom_boxplot(aes(fill = Category)) +  # if you have "Source"
  facet_grid(. ~ Category) +            # vertical split by source
  theme_minimal() +
  labs(
    title = "ARG per 16S rRNA gene ratio per Antimicrobial Class",
    x = "ARGs per 16S copy",
    y = "Gene (Group)"
  )  +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),  # facet labels
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

##both times together

# Add time assignment (you'll need to adjust these sample lists based on your study design)
arg_summary <- arg_summary %>%
  mutate(
    Time = case_when(
      # Define Time 1 samples (adjust these sample IDs based on your data)
      Sample %in% c("E06", "E02", "E10", "F01", "F04","E03", "E07", "F05", "F09", "E01", "E05", "E09", "E12", "F03", "F11", "E04", "E08", "E11", "F02", "F06", "F10" ) ~ "Time 1",
      # Define Time 2 samples (adjust these sample IDs based on your data)  
      Sample %in% c("F12", "G04", "G09", "G12", "H04", "H08","H01", "H05", "H09",  "G03", "G08", "G11", "H03", "H07", "G02", "G06", "G10", "H06", "H10") ~ "Time 2",
      # Add more sample assignments as needed
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(Time != "Unknown")  # Remove samples without time assignment


# Create a combined variable for better plotting
arg_summary <- arg_summary %>%
  mutate(
    Category_Time = paste(Category, "-", Time),
    Time_Factor = factor(Time, levels = c("Time 1", "Time 2"))
  )

# Create the plot with time as color and faceted by category
ggplot(arg_summary, aes(x = ARG_per_16S, y = Class, fill = Time_Factor)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Time 1" = "blue", "Time 2" = "red")) +
  facet_grid(. ~ Category) +
  theme_minimal() +
  labs(
    title = "ARG per 16S rRNA gene ratio per Antimicrobial Class by Time Point",
    x = "ARG per copy of 16S rRNA",
    y = "Gene class",
    fill = "Time Point"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "bottom"
  )

##deleting the underscore 
ggplot(arg_summary, aes(x = ARG_per_16S, y = Class, fill = Time_Factor)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Time 1" = "darkblue", "Time 2" = "darkred")) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +  # Add this line
  facet_grid(. ~ Category) +
  theme_minimal() +
  labs(
   # title = "ARG per 16S rRNA gene ratio per Antimicrobial Class by Time Point",
    x = "ARG per copy of 16S rRNA",
    y = "Gene class",
    fill = "Time Point"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "right"
  )
## same panel

# Plot with time as facets and category as fill
ggplot(arg_summary, aes(x = ARG_per_16S, y = Class, fill = Time_Factor)) +
  scale_fill_manual(values = c("Time 1" = "darkblue", "Time 2" = "darkred")) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +  # Add this line
  geom_boxplot(alpha = 0.7) +
  facet_grid(Time_Factor ~ Category) +
  theme_minimal() +
  labs(
  #  title = "ARG per 16S rRNA gene ratio per Antimicrobial Class by Time Point",
    x = "ARG per copy of 16S rRNA",
    y = "Antimicorbial class",
    fill = "Sample time"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )
