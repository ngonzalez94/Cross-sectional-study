# =============================================================================
# Resistome PCoA Analysis
# Classes: Beta-lactams, Tetracyclines, MLS Multidrug Resistance, Aminoglycosides
# Distance: Bray-Curtis + Jaccard
# Shapes: Time 1 vs Time 2 | Colours: 4 sources
# =============================================================================

# --- 1. INSTALL & LOAD PACKAGES ----------------------------------------------

# Uncomment to install if needed:
# install.packages(c("tidyverse", "vegan", "ggplot2", "patchwork"))

library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)

# --- 2. LOAD YOUR DATA -------------------------------------------------------
# Adjust file paths to where your files are located
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA/")
milk_data   <- read_csv("milk_compiled.csv")    # <-- change filename
faeces_data <- read_csv("faeces_compiled.csv")  # <-- change filename

# Combine both files into one dataframe
combined <- bind_rows(milk_data, faeces_data)

# Expected columns: Sample, Gene, Total_Hits
# e.g.: E03  MEG_1189|Drugs|Bacitracin|...|BACA  45

# --- 3. PARSE GENE ANNOTATION ------------------------------------------------
# Gene format: MEG_ID|Drugs|Class|Subclass|GeneName
# We extract: Class (position 3) and GeneName (position 5)

combined <- combined %>%
  mutate( 
    Level = str_split_fixed(Gene, "\\|", 5)[, 2]
   )

# --- 4. SAMPLE METADATA ------------------------------------------------------

# 4a. Time point assignment
# Naming: E01-F11 = Time 1 | F12-H10 = Time 2
# Logic: extract letter prefix and numeric suffix, then classify

get_timepoint <- function(sample_id) {
  letter <- str_extract(sample_id, "^[A-Z]+")
  number <- as.integer(str_extract(sample_id, "[0-9]+$"))
  
  case_when(
    letter == "E"                          ~ "Time 1",
    letter == "F" & number <= 11           ~ "Time 1",
    letter == "F" & number >= 12           ~ "Time 2",
    letter %in% c("G", "H")               ~ "Time 2",
    TRUE                                   ~ NA_character_
  )
}

# 4b. Source assignment
get_source <- function(sample_id) {
  case_when(
    sample_id %in% c("E06","E02","E10","F01","F04","F12",
                     "G04","G09","G12","H04","H08")        ~ "Keeper Faeces",
    sample_id %in% c("E03","E07","F05","F09","G01",
                     "G05","H01","H05","H09")              ~ "Replacement Milk",
    sample_id %in% c("E01","E05","E09","E12","F03","F11",
                     "G03","G08","G11","H03","H07")        ~ "Bobby Faeces",
    sample_id %in% c("E04","E08","E11","F02","F06","F10",
                     "G02","G06","G10","H06","H10")        ~ "Waste Milk",
    TRUE                                                   ~ NA_character_
  )
}

# Build metadata table (one row per unique sample)
metadata <- combined %>%
  distinct(Sample) %>%
  mutate(
    TimePoint = get_timepoint(Sample),
    Source    = get_source(Sample)
  )
# --- 5. DEFINE ANTIBIOTIC CLASSES OF INTEREST --------------------------------

# Check what class names exist in your data:
 unique(combined$Level)
# Quick check — print any samples not assigned to a source or timepoint
unassigned <- metadata %>% filter(is.na(TimePoint) | is.na(Source))
if (nrow(unassigned) > 0) {
  message("⚠️  Unassigned samples found:")
  print(unassigned)
} else {
  message("✅  All samples assigned to TimePoint and Source.")
}

# Map your data's class labels to the 4 classes of interest
# Adjust the strings below if they differ in your data!

classes_of_interest <- c(
  "Drugs"
)


# --- 6. PLOTTING SETTINGS ----------------------------------------------------

source_colours <- c(
  "Keeper Faeces"    = "#377EB8",   
  "Replacement Milk" = "#4DAF4A",   
  "Bobby Faeces"     = "#E41A1C",   
  "Waste Milk"       = "#984EA3"    
)

time_shapes <- c(
  "Time 1" = 16,   # filled circle
  "Time 2" = 17    # filled triangle
)

# --- 7. CORE PCoA FUNCTION ---------------------------------------------------

run_pcoa <- function(data, meta, ab_class, distance_method) {
  # Filter to the class of interest
  class_data <- data %>%
    filter(Level == ab_class) %>%
    select(Sample, GeneName, Total_Hits)
  
  # Aggregate in case of duplicate Sample+Gene rows
  class_data <- class_data %>%
    group_by(Sample, GeneName) %>%
    summarise(Total_Hits = sum(Total_Hits), .groups = "drop")
  
  # Pivot to wide matrix: rows = samples, cols = genes
  mat_wide <- class_data %>%
    pivot_wider(names_from = GeneName, values_from = Total_Hits, values_fill = 0) %>%
    column_to_rownames("Sample")
  
  # Keep only samples present in metadata
  shared_samples <- intersect(rownames(mat_wide), meta$Sample)
  
  if (length(shared_samples) < 3) {
    message("⚠️  Not enough samples for class: ", ab_class, " (n=", length(shared_samples), ")")
    return(NULL)
  }
  
  mat_wide <- mat_wide[shared_samples, , drop = FALSE]
  
  # For Jaccard: convert to presence/absence
  if (distance_method == "jaccard") {
    mat_input <- decostand(mat_wide, method = "pa")  # presence/absence
  } else {
    mat_input <- mat_wide  # raw counts for Bray-Curtis
  }
  
  # Remove all-zero genes (columns) after potential PA conversion
  mat_input <- mat_input[, colSums(mat_input) > 0, drop = FALSE]
  
  # Calculate distance matrix
  dist_mat <- vegdist(mat_input, method = distance_method, binary = (distance_method == "jaccard"))
  
  # Run PCoA
  pcoa_result <- cmdscale(dist_mat, k = 2, eig = TRUE)
  
  # Extract % variance explained
  eig_vals   <- pcoa_result$eig
  eig_pos    <- eig_vals[eig_vals > 0]
  var_exp    <- round(eig_pos[1:2] / sum(eig_pos) * 100, 1)
  
  # Build plot dataframe
  plot_df <- as.data.frame(pcoa_result$points) %>%
    rownames_to_column("Sample") %>%
    rename(PC1 = V1, PC2 = V2) %>%
    left_join(meta, by = "Sample")
  
  # Build plot
  dist_label <- ifelse(distance_method == "bray", "Bray-Curtis", "Jaccard")
  
  p <- ggplot(plot_df, aes(x = PC1, y = PC2,
                           colour = Source,
                           shape  = TimePoint)) +
    geom_point(size = 4, alpha = 0.85) +
    stat_ellipse(aes(group = Source, colour = Source),
                 type = "t", level = 0.75,
                 linetype = "dashed", linewidth = 0.5,
                 show.legend = FALSE) +
    scale_colour_manual(values = source_colours, 
                        labels = c("Bobby Faeces" = "CMR-fed calves' faeces", "Keeper Faeces" = "WM-fed calves' faeces"))+
    scale_shape_manual(values = time_shapes) +
    labs(
      #title    = paste0(ab_class, " — ", dist_label),
      x        = paste0("PC1 (", var_exp[1], "%)"),
      y        = paste0("PC2 (", var_exp[2], "%)"),
      colour   = "Source",
      shape    = "Time Point"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# --- 8. GENERATE ALL PLOTS ---------------------------------------------------

all_plots <- list()

for (ab_class in classes_of_interest) {
  p_bc  <- run_pcoa(combined, metadata, ab_class, "bray")
  p_jac <- run_pcoa(combined, metadata, ab_class, "jaccard")
  
  if (!is.null(p_bc) && !is.null(p_jac)) {
    # Side-by-side: Bray-Curtis | Jaccard
    combined_plot <- p_bc + p_jac +
      plot_layout(guides = "collect") &
      theme(legend.position = "right")
    
    all_plots[[ab_class]] <- combined_plot
  }
}

# --- 9. DISPLAY PLOTS --------------------------------------------------------

# Option A: Print all plots to your R graphics device (e.g. RStudio viewer)
for (ab_class in names(all_plots)) {
  print(all_plots[[ab_class]])
}

# Option B: Save each class as its own PDF page (one file per class)
for (ab_class in names(all_plots)) {
  safe_name <- str_replace_all(ab_class, "[^A-Za-z0-9]", "_")
  ggsave(
    filename = paste0("PCoA_", safe_name, ".pdf"),
    plot     = all_plots[[ab_class]],
    width    = 14,
    height   = 6
  )
  message("Saved: PCoA_", safe_name, ".pdf")
}

# Option C: Save ALL classes together in one multi-page PDF
pdf("PCoA_all_classes.pdf", width = 14, height = 6)
for (ab_class in names(all_plots)) {
  print(all_plots[[ab_class]])
}
dev.off()
message("✅  Saved: PCoA_all_classes.pdf")

# --- 10. OPTIONAL: PERMANOVA (adonis2) ----------------------------------------
# Tests whether source / timepoint significantly explain community composition
# Uncomment to run:

for (ab_class in classes_of_interest) {
  
  class_data <- combined %>%
    filter(Level == ab_class) %>%
    group_by(Sample, GeneName) %>%
    summarise(Total_Hits = sum(Total_Hits), .groups = "drop") %>%
    pivot_wider(names_from  = GeneName,
                values_from = Total_Hits,
                values_fill = 0) %>%
    column_to_rownames("Sample")
  
  # Match samples present in BOTH matrix and metadata
  shared <- intersect(rownames(class_data), metadata$Sample)
  
  if (length(shared) < 3) {
    message("Skipping ", ab_class, ": too few samples (n=", length(shared), ")")
    next
  }
  
  class_data <- class_data[shared, , drop = FALSE]
  
  # ✅ Key fix: filter AND reorder metadata to match matrix row order exactly
  meta_sub <- metadata %>%
    filter(Sample %in% shared) %>%
    arrange(match(Sample, shared))   # <-- this is what was missing
  
  dist_bc <- vegdist(class_data, method = "bray")
  
  perm <- adonis2(dist_bc ~ Source + TimePoint,
                  data         = meta_sub,
                  permutations = 999)
  
  cat("\n=== PERMANOVA:", ab_class, "===\n")
  print(perm)
}


# if your overall PERMANOVA is significant for Source, you'll want to know which sources differ from each other:
# install.packages("pairwiseAdonis")  # from GitHub
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis2(dist_bc ~ Source, data = meta_sub, permutations = 999)
pairwise.adonis2(dist_bc ~ TimePoint, data = meta_sub, permutations = 999)
