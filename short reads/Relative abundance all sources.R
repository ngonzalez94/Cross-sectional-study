# Load required libraries
library(tidyverse)
install.packages("tibble")
library(tibble)
library(readxl)
library(dplyr)

# Define the file path
file_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken/Relative abundance.xlsx"

# Load the bracken data and convert to tibble
dbf <- read_excel(file_path, sheet=2) %>% as_tibble()
dkf <- read_excel(file_path, sheet=4) %>% as_tibble()
drm <- read_excel(file_path, sheet=8) %>% as_tibble()
dwm <- read_excel(file_path, sheet=6) %>% as_tibble()

# Function to prepare data and get top taxa
prepare_data <- function(df, sample_name) {
  df <- as.data.frame(df)
  taxon_col <- colnames(df)[1]
  
  t1_cols <- grep("_T1$", colnames(df), value = TRUE)
  t2_cols <- grep("_T2$", colnames(df), value = TRUE)
  
  df_summary <- df %>%
    as_tibble() %>%
    mutate(
      Time1 = rowMeans(select(., all_of(t1_cols)), na.rm = TRUE),
      Time2 = rowMeans(select(., all_of(t2_cols)), na.rm = TRUE)
    ) %>%
    select(taxon = 1, Time1, Time2) %>%
    filter(!is.na(Time1) | !is.na(Time2)) %>%
    mutate(
      Time1 = replace_na(Time1, 0),
      Time2 = replace_na(Time2, 0)
    )
  
  total_t1 <- sum(df_summary$Time1, na.rm = TRUE)
  total_t2 <- sum(df_summary$Time2, na.rm = TRUE)
  
  df_summary <- df_summary %>%
    mutate(
      Time1 = Time1 / total_t1,
      Time2 = Time2 / total_t2
    )
  
  # Process Time1
  top10_t1 <- df_summary %>%
    filter(Time1 > 0) %>%
    arrange(desc(Time1)) %>%
    slice(1:10) %>%
    pull(taxon)
  
  df_t1 <- df_summary %>%
    mutate(taxon_group = ifelse(taxon %in% top10_t1, as.character(taxon), "Other")) %>%
    group_by(taxon_group) %>%
    summarise(Abundance = sum(Time1, na.rm = TRUE), .groups = "drop") %>%
    rename(taxon = taxon_group) %>%
    mutate(Time = "Time1", Sample = sample_name)
  
  # Process Time2
  top10_t2 <- df_summary %>%
    filter(Time2 > 0) %>%
    arrange(desc(Time2)) %>%
    slice(1:10) %>%
    pull(taxon)
  
  df_t2 <- df_summary %>%
    mutate(taxon_group = ifelse(taxon %in% top10_t2, as.character(taxon), "Other")) %>%
    group_by(taxon_group) %>%
    summarise(Abundance = sum(Time2, na.rm = TRUE), .groups = "drop") %>%
    rename(taxon = taxon_group) %>%
    mutate(Time = "Time2", Sample = sample_name)
  
  df_final <- bind_rows(df_t1, df_t2) %>%
    mutate(Time = factor(Time, levels = c("Time1", "Time2"))) %>%
    filter(!is.na(Abundance), Abundance > 0)
  
  return(df_final)
}

# Combine all data
combined_data <- bind_rows(
  prepare_data(dbf, "BF"),
  prepare_data(dkf, "KF"),
  prepare_data(drm, "RM"),
  prepare_data(dwm, "WM")
)

# Order samples
combined_data$Sample <- factor(combined_data$Sample, 
                               levels = c("BF", "KF", "RM", "WM"))

# Calculate overall abundance for each taxon to order legend
taxon_order <- combined_data %>%
  group_by(taxon) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  filter(taxon != "Other") %>%
  pull(taxon)

# Add "Other" at the end
taxon_order <- c(taxon_order, "Other")

# Set factor levels for taxa (controls legend order)
combined_data$taxon <- factor(combined_data$taxon, levels = rev(taxon_order))

# DEFINE YOUR CUSTOM COLORS HERE (24 colors for taxa + grey for "Other")
custom_colors <- c(
"#d9ef8b",
"#a6d96a",
"#66bd63",
"#1a9850",
"#006890",
"#40004b",
"#762a83",
"#9970ab",
"#c2a5cf",
"#e7d4e8",
"#5aae61",
"#1b7837",
"#00441b",
"#a50026",
"#d73027",
"#f46d43",
"#fdae61",
"#fee090",
"#ffffbf",
"#e0f3f8",
"#abd9e9",
"#74add1",
"#4575b4",
"#313695"
)

# Get list of taxa (excluding "Other")
taxa_list <- levels(combined_data$taxon)
taxa_no_other <- setdiff(taxa_list, "Other")

# Create color vector
color_vector <- setNames(
  c("grey70", custom_colors[1:length(taxa_no_other)]),
  c("Other", taxa_no_other)
)

# Create the plot
ggplot(combined_data, aes(x = Sample, y = Abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_vector) +
  facet_wrap(~Time, ncol = 1, strip.position = "right", 
             labeller = labeller(Time = c(Time1 = "Time 1", Time2 = "Time 2"))) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(labels = c(
    "BF" = "CMR_f",
    "KF" = "WM_f",
    "RM" = "CMR",
    "Wk" = "WM"
  )) +
  theme_minimal() +
  labs(
    x = "Sample source",
    y = "Relative Abundance",
    fill = "Family"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 16)
  )
##STATISTICS##

# -----------------------------------------------------------
# 1. Reshape all sheets into one long dataframe
#    extracting farm, sample type, and timepoint from column names
# -----------------------------------------------------------
reshape_sheet <- function(df, sample_type) {
  df %>%
    rename(taxon = 1) %>%
    pivot_longer(cols = -taxon, names_to = "sample_id", values_to = "abundance") %>%
    mutate(
      farm        = str_extract(sample_id, "F[0-9]+"),
      time        = str_extract(sample_id, "T[12]$"),
      sample_type = sample_type
    ) %>%
    filter(!is.na(abundance))
}

raw_long <- bind_rows(
  reshape_sheet(dbf, "Bobby"),
  reshape_sheet(dkf, "Keeper"),
  reshape_sheet(drm, "RawMilk"),
  reshape_sheet(dwm, "WholeMilk")
) %>%
  # Normalise to relative abundance within each sample
  group_by(sample_type, sample_id) %>%
  mutate(rel_abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    sample_type = factor(sample_type, levels = c("Bobby", "Keeper", "RawMilk", "WholeMilk")),
    time        = factor(time, levels = c("T1", "T2")),
    farm        = factor(farm)
  )

# -----------------------------------------------------------
# Get the plotted taxa from your existing combined_data
# (run your original plot script first so combined_data exists)
# -----------------------------------------------------------
plotted_taxa <- combined_data %>%
  filter(taxon != "Other") %>%
  pull(taxon) %>%
  as.character() %>%
  unique()

raw_long <- raw_long %>% filter(taxon %in% plotted_taxa)

# -----------------------------------------------------------
# STAGE 1: Time effect
# Within each sample type, does T1 differ from T2?
# Wilcoxon signed-rank paired by farm (Farm1_T1 vs Farm1_T2)
# -----------------------------------------------------------
stage1_results <- expand_grid(
  taxon       = plotted_taxa,
  sample_type = levels(raw_long$sample_type)
) %>%
  mutate(
    wx = map2(taxon, sample_type, function(tx, st) {
      paired_df <- raw_long %>%
        filter(taxon == tx, sample_type == st) %>%
        select(farm, time, rel_abundance) %>%
        pivot_wider(id_cols = farm, names_from = time, values_from = rel_abundance)
      
      if (!all(c("T1", "T2") %in% colnames(paired_df))) return(NULL)
      paired_df <- paired_df %>% filter(!is.na(T1), !is.na(T2))
      if (nrow(paired_df) < 2) return(NULL)
      
      tryCatch({
        wt <- wilcox.test(paired_df$T1, paired_df$T2,
                          paired = TRUE, exact = FALSE, correct = TRUE)
        z_stat  <- qnorm(wt$p.value / 2) * sign(median(paired_df$T1 - paired_df$T2))
        n_pairs <- nrow(paired_df)
        tibble(
          n_pairs   = n_pairs,
          W         = wt$statistic,
          Z         = z_stat,
          p_raw     = wt$p.value,
          r_effect  = abs(z_stat) / sqrt(n_pairs),
          direction = ifelse(median(paired_df$T1) > median(paired_df$T2), "T1 > T2", "T2 > T1")
        )
      }, error = function(e) NULL)
    })
  ) %>%
  filter(!map_lgl(wx, is.null)) %>%
  unnest(wx) %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH"))

cat("\n===== STAGE 1: Time effect (T1 vs T2 within each sample type) =====\n")
# the test itself asks "across all 6 farms, is there a consistent T1 vs T2 difference in Bobby samples?"
stage1_results %>%
  arrange(p_adj) %>%
  filter(p_adj < 0.05) %>%
  print(n = 50)

# Which sample types have NO time effect (safe to collapse)?
no_time_effect <- stage1_results %>%
  group_by(sample_type) %>%
  summarise(any_sig = any(p_adj < 0.05), .groups = "drop") %>%
  filter(!any_sig) %>%
  pull(sample_type)

cat("\nSample types with no time effect (T1+T2 will be collapsed):\n")
print(no_time_effect)

time_effect <- setdiff(levels(raw_long$sample_type), no_time_effect)
cat("Sample types WITH time effect (kept separate):\n")
print(time_effect)

# -----------------------------------------------------------
# STAGE 2: Farm effect
# Within each sample type, do farms differ from each other?
# Collapse T1+T2 for types with no time effect
# Kruskal-Wallis across farms, Dunn post-hoc if significant
# -----------------------------------------------------------

# Collapse timepoints for sample types with no time effect
raw_collapsed <- raw_long %>%
  filter(sample_type %in% no_time_effect) %>%
  group_by(taxon, sample_type, farm) %>%
  summarise(rel_abundance = mean(rel_abundance, na.rm = TRUE), .groups = "drop")

# Keep separate timepoints for sample types with time effect
raw_stratified <- raw_long %>%
  filter(sample_type %in% time_effect)

stage2_results <- bind_rows(
  # Collapsed (no time effect)
  expand_grid(
    taxon       = plotted_taxa,
    sample_type = no_time_effect
  ) %>%
    mutate(time = "collapsed"),
  
  # Stratified (time effect present)
  expand_grid(
    taxon       = plotted_taxa,
    sample_type = time_effect,
    time        = c("T1", "T2")
  )
) %>%
  mutate(
    kw = pmap(list(taxon, sample_type, time), function(tx, st, tm) {
      if (tm == "collapsed") {
        sub_df <- raw_collapsed %>% filter(taxon == tx, sample_type == st)
      } else {
        sub_df <- raw_stratified %>% filter(taxon == tx, sample_type == st, time == tm)
      }
      if (length(unique(sub_df$farm)) < 2) return(NULL)
      tryCatch(
        kruskal.test(rel_abundance ~ farm, data = sub_df),
        error = function(e) NULL
      )
    })
  ) %>%
  filter(!map_lgl(kw, is.null)) %>%
  mutate(
    kw_statistic = map_dbl(kw, ~ .x$statistic),
    kw_df        = map_dbl(kw, ~ .x$parameter),
    kw_p         = map_dbl(kw, ~ .x$p.value)
  ) %>%
  select(-kw) %>%
  mutate(kw_p_adj = p.adjust(kw_p, method = "BH"))

cat("\n===== STAGE 2: Farm effect (within each sample type) =====\n")
stage2_results %>%
  filter(kw_p_adj < 0.05) %>%
  arrange(sample_type, kw_p_adj) %>%
  print(n = 50)

# Dunn post-hoc for significant farm effects
stage2_dunn <- stage2_results %>%
  filter(kw_p_adj < 0.05) %>%
  select(taxon, sample_type, time) %>%
  mutate(
    dunn = pmap(list(taxon, sample_type, time), function(tx, st, tm) {
      if (tm == "collapsed") {
        sub_df <- raw_collapsed %>% filter(taxon == tx, sample_type == st)
      } else {
        sub_df <- raw_stratified %>% filter(taxon == tx, sample_type == st, time == tm)
      }
      tryCatch({
        dt <- dunn.test::dunn.test(
          x = sub_df$rel_abundance, g = sub_df$farm,
          method = "bh", kw = FALSE, label = TRUE, altp = TRUE
        )
        tibble(
          comparison = dt$comparisons,
          Z          = dt$Z,
          p_raw      = dt$altP,
          p_adj      = dt$altP.adjusted
        )
      }, error = function(e) NULL)
    })
  ) %>%
  filter(!map_lgl(dunn, is.null)) %>%
  unnest(dunn)

cat("\n--- Stage 2 Dunn post-hoc: which farms differ? ---\n")
stage2_dunn %>%
  filter(p_adj < 0.05) %>%
  arrange(sample_type, taxon, p_adj) %>%
  print(n = 100)

# Which sample types have NO farm effect?
no_farm_effect <- stage2_results %>%
  group_by(sample_type) %>%
  summarise(any_sig = any(kw_p_adj < 0.05), .groups = "drop") %>%
  filter(!any_sig) %>%
  pull(sample_type)

cat("\nSample types with no farm effect (farms will be collapsed for stage 3):\n")
print(no_farm_effect)

# -----------------------------------------------------------
# STAGE 3: Group effect
# Only for sample types with no farm effect
# Compare Bobby vs Keeper (calf groups)
# Compare RawMilk vs WholeMilk (milk groups)
# Wilcoxon rank-sum (unpaired, independent groups)
# -----------------------------------------------------------

# Collapse farms for sample types with no farm effect
raw_group <- raw_long %>%
  filter(sample_type %in% no_farm_effect) %>%
  {
    if (length(no_time_effect) > 0)
      group_by(., taxon, sample_type, farm) %>%
      summarise(rel_abundance = mean(rel_abundance, na.rm = TRUE), .groups = "drop")
    else .
  }

# Calf comparison: Bobby vs Keeper
calf_types <- intersect(no_farm_effect, c("Bobby", "Keeper"))
milk_types  <- intersect(no_farm_effect, c("RawMilk", "WholeMilk"))

run_group_comparison <- function(taxa, type_a, type_b, data) {
  map_dfr(taxa, function(tx) {
    grp_a <- data %>% filter(taxon == tx, sample_type == type_a) %>% pull(rel_abundance)
    grp_b <- data %>% filter(taxon == tx, sample_type == type_b) %>% pull(rel_abundance)
    if (length(grp_a) < 2 || length(grp_b) < 2) return(NULL)
    tryCatch({
      wt <- wilcox.test(grp_a, grp_b, paired = FALSE, exact = FALSE, correct = TRUE)
      z_stat <- qnorm(wt$p.value / 2) * sign(median(grp_a) - median(grp_b))
      tibble(
        taxon      = tx,
        comparison = paste(type_a, "vs", type_b),
        W          = wt$statistic,
        Z          = z_stat,
        p_raw      = wt$p.value,
        r_effect   = abs(z_stat) / sqrt(length(grp_a) + length(grp_b)),
        direction  = ifelse(median(grp_a) > median(grp_b),
                            paste(type_a, ">", type_b),
                            paste(type_b, ">", type_a))
      )
    }, error = function(e) NULL)
  })
}

stage3_results <- bind_rows(
  if (length(calf_types) == 2)
    run_group_comparison(plotted_taxa, "Bobby", "Keeper", raw_group),
  if (length(milk_types) == 2)
    run_group_comparison(plotted_taxa, "RawMilk", "WholeMilk", raw_group)
) %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH"))

cat("\n===== STAGE 3: Group effect =====\n")
stage3_results %>%
  filter(p_adj < 0.05) %>%
  arrange(comparison, p_adj) %>%
  print(n = 50)

# -----------------------------------------------------------
# Export all results
# -----------------------------------------------------------
output_dir <- dirname(file_path)
write_csv(stage1_results, file.path(output_dir, "stage1_time_effect.csv"))
write_csv(stage2_results, file.path(output_dir, "stage2_farm_effect.csv"))
write_csv(stage2_dunn,    file.path(output_dir, "stage2_farm_dunn.csv"))
write_csv(stage3_results, file.path(output_dir, "stage3_group_effect.csv"))
cat("\nAll results saved to:", output_dir, "\n")


