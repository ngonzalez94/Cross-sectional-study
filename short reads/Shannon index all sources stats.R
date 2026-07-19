library(vegan)
library(tidyverse)

# ============================================================
# Calculate Alpha Diversity Metrics
# ============================================================

# Define the file path
file_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken/Relative abundance.xlsx"

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
rm_data <- prepare_pcoa_data(drm, "R")
wm_data <- prepare_pcoa_data(dwm, "W")

# Combine all data
all_data <- bind_rows(bf_data, kf_data, rm_data, wm_data)

# Separate metadata from abundance data
metadata <- all_data %>%
  select(Sample, Time, Group)

abundance_data <- all_data %>%
  select(-Sample, -Time, -Group)

# Replace NAs with 0
abundance_data[is.na(abundance_data)] <- 0

# Calculate multiple alpha diversity indices
alpha_diversity <- data.frame(
  Sample = metadata$Sample,
  Group = metadata$Group,
  Time = metadata$Time,
  
  # Species richness (number of taxa present)
  Richness = apply(abundance_data > 0, 1, sum),
  
  # Shannon diversity index (accounts for abundance and evenness)
  Shannon = diversity(abundance_data, index = "shannon"),
  
  # Simpson diversity index (probability two individuals are different species)
  Simpson = diversity(abundance_data, index = "simpson"),
  
  # Inverse Simpson (more intuitive - higher = more diverse)
  InvSimpson = diversity(abundance_data, index = "invsimpson"),
  
  # Pielou's evenness (Shannon / log(Richness))
  Pielou = diversity(abundance_data, index = "shannon") / log(apply(abundance_data > 0, 1, sum))
)

# View the results
print(alpha_diversity)

# Summary statistics by group
cat("\n=== Alpha Diversity by Sample Type ===\n")
alpha_diversity %>%
  group_by(Group) %>%
  summarise(
    Mean_Richness = mean(Richness),
    SD_Richness = sd(Richness),
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Simpson = mean(Simpson),
    SD_Simpson = sd(Simpson)
  ) %>%
  print()

cat("\n=== Alpha Diversity by Time Point ===\n")
alpha_diversity %>%
  group_by(Time) %>%
  summarise(
    Mean_Richness = mean(Richness),
    SD_Richness = sd(Richness),
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Simpson = mean(Simpson),
    SD_Simpson = sd(Simpson)
  ) %>%
  print()

# ============================================================
# Visualize Alpha Diversity
# ============================================================

# 1. Shannon Diversity Boxplot
p_shannon <- ggplot(alpha_diversity, aes(x = Group, y = Shannon, fill = Time)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.75)) +
  geom_point(aes(shape = Time), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75),
                                             size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Time1" = "darkblue", "Time2" = "darkred"),
                    labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_shape_manual(values = c("Time1" = 16, "Time2" = 16),
                     labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_x_discrete(labels = c("B" = "CMR_f", "K" = "WM_f", "RM" = "CMR", "WM" = "WM")) +
  theme_minimal() +
  labs(
    x = "Sample Type",
    y = "Shannon Index",
    fill = "Time Point",
    shape = "Time Point"
  ) +
  theme(legend.position = "right")

# 2. Species Richness Boxplot
p_richness <- ggplot(alpha_diversity, aes(x = Group, y = Richness, fill = Time)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.75)) +
  geom_point(aes(shape = Time), 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.75),
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Time1" = "darkblue", "Time2" = "darkred"),
                    labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_shape_manual(values = c("Time1" = 16, "Time2" = 16),
                     labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_x_discrete(labels = c("B" = "CMR_f", "K" = "WM_f", "RM" = "CMR", "WM" = "WM")) +
  theme_minimal() +
  labs(
    #title = "Chao1",
    x = "Sample Type",
    y = "Chao 1 diversity",
    fill = "Time Point",
    shape = "Time Point"
  ) +
  theme(legend.position = "right")

# 3. Simpson Diversity Boxplot
p_simpson <- ggplot(alpha_diversity, aes(x = Group, y = Simpson, fill = Time)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.75)) +
  geom_point(aes(shape = Time), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75),
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Time1" = "darkblue", "Time2" = "darkred"),
                    labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_shape_manual(values = c("Time1" = 16, "Time2" = 16),
                     labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_x_discrete(labels = c("B" = "CMR_f", "K" = "WM_f", "RM" = "CMR", "WM" = "WM")) +
  theme_minimal() +
  labs(
    x = "Sample Type",
    y = "Simpson Index",
    fill = "Time Point",
    shape = "Time Point"
  ) +
  theme(legend.position = "right")

# Display plots
print(p_shannon)
print(p_richness)
print(p_simpson)

# Or combine all three
library(patchwork)
(p_shannon / p_simpson /p_richness)+
  plot_layout(guides ="collect")

# ============================================================
# Statistical Tests
# ============================================================

# Test for differences between groups using ANOVA or Kruskal-Wallis

# Shannon diversity
cat("\n=== Statistical Tests for Shannon Diversity ===\n")

# Test for normality
shapiro_shannon <- shapiro.test(alpha_diversity$Shannon)
cat("Shapiro-Wilk test for normality (p > 0.05 = normal):", shapiro_shannon$p.value, "\n")

if(shapiro_shannon$p.value > 0.05) {
  # Use ANOVA if normal
  aov_group <- aov(Shannon ~ Group, data = alpha_diversity)
  cat("\nANOVA - Effect of Sample Type:\n")
  print(summary(aov_group))
  
  # Post-hoc test
  cat("\nTukey HSD post-hoc test:\n")
  print(TukeyHSD(aov_group))
  
  # Test time effect
  aov_time <- aov(Shannon ~ Time, data = alpha_diversity)
  cat("\nANOVA - Effect of Time:\n")
  print(summary(aov_time))
  
  # Test interaction
  aov_both <- aov(Shannon ~ Group * Time, data = alpha_diversity)
  cat("\nANOVA - Effect of Sample Type and Time:\n")
  print(summary(aov_both))
  
} else {
  # Use Kruskal-Wallis if not normal
  kw_group <- kruskal.test(Shannon ~ Group, data = alpha_diversity)
  cat("\nKruskal-Wallis test - Effect of Sample Type:\n")
  print(kw_group)
  
  # Post-hoc test (pairwise Wilcoxon)
  cat("\nPairwise Wilcoxon tests:\n")
  print(pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Group, p.adjust.method = "BH"))
  
  # Test time effect
  kw_time <- kruskal.test(Shannon ~ Time, data = alpha_diversity)
  cat("\nKruskal-Wallis test - Effect of Time:\n")
  print(kw_time)
}

# Similar tests for Richness
cat("\n\n=== Statistical Tests for Species Richness ===\n")
shapiro_richness <- shapiro.test(alpha_diversity$Richness)
cat("Shapiro-Wilk test for normality:", shapiro_richness$p.value, "\n")

if(shapiro_richness$p.value > 0.05) {
  aov_richness <- aov(Richness ~ Group * Time, data = alpha_diversity)
  print(summary(aov_richness))
} else {
  kw_richness <- kruskal.test(Richness ~ Group, data = alpha_diversity)
  print(kw_richness)
}

# ============================================================
# Hutcheson t-test for Shannon Diversity
# ============================================================

# Function to perform Hutcheson t-test
hutcheson_t_test <- function(H1, H2, N1, N2, S1, S2) {
  # H1, H2 = Shannon diversity indices
  # N1, N2 = Total number of individuals (sum of abundances)
  # S1, S2 = Number of species (richness)
  
  # Calculate variance of Shannon index for each sample
  var_H1 <- (sum(sapply(1:S1, function(i) {
    pi <- abundance_data[sample1_idx, i] / N1
    if(pi > 0) pi * (log(pi))^2 else 0
  })) - H1^2) / N1 + (S1 - 1) / (2 * N1^2)
  
  var_H2 <- (sum(sapply(1:S2, function(i) {
    pi <- abundance_data[sample2_idx, i] / N2
    if(pi > 0) pi * (log(pi))^2 else 0
  })) - H2^2) / N2 + (S2 - 1) / (2 * N2^2)
  
  # Calculate t-statistic
  t_stat <- abs(H1 - H2) / sqrt(var_H1 + var_H2)
  
  # Degrees of freedom
  df <- (var_H1 + var_H2)^2 / ((var_H1^2 / N1) + (var_H2^2 / N2))
  
  # P-value (two-tailed)
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  
  return(list(t_statistic = t_stat, df = df, p_value = p_value))
}

# However, there's a simpler approach using the BiodiversityR package
# or we can use a custom implementation

# ============================================================
# Method 1: Using a custom Hutcheson test function
# ============================================================

hutcheson_test <- function(data, group_var, group1, group2) {
  # Subset data for two groups
  data1 <- data[data[[group_var]] == group1, ]
  data2 <- data[data[[group_var]] == group2, ]
  
  # Calculate mean Shannon and sample sizes
  H1 <- mean(data1$Shannon)
  H2 <- mean(data2$Shannon)
  N1 <- nrow(data1)
  N2 <- nrow(data2)
  
  # Calculate variance of Shannon
  var_H1 <- var(data1$Shannon)
  var_H2 <- var(data2$Shannon)
  
  # Calculate t-statistic
  t_stat <- (H1 - H2) / sqrt((var_H1/N1) + (var_H2/N2))
  
  # Degrees of freedom (Welch-Satterthwaite equation)
  df <- ((var_H1/N1) + (var_H2/N2))^2 / 
    ((var_H1/N1)^2/(N1-1) + (var_H2/N2)^2/(N2-1))
  
  # P-value (two-tailed)
  p_value <- 2 * pt(-abs(t_stat), df)
  
  result <- data.frame(
    Group1 = group1,
    Group2 = group2,
    Mean_H1 = H1,
    Mean_H2 = H2,
    Difference = H1 - H2,
    t_statistic = t_stat,
    df = df,
    p_value = p_value,
    Significant = ifelse(p_value < 0.05, "Yes", "No")
  )
  
  return(result)
}

# ============================================================
# Compare Shannon Diversity Between Groups
# ============================================================

cat("=== Hutcheson t-test for Shannon Diversity ===\n\n")

# 1. Compare between Time points
cat("--- Time 1 vs Time 2 ---\n")
hutcheson_time <- hutcheson_test(alpha_diversity, "Time", "Time1", "Time2")
print(hutcheson_time)

# 2. Compare between Sample Types (pairwise)
cat("\n--- Pairwise Comparisons Between Sample Types ---\n")

groups <- unique(alpha_diversity$Group)
comparisons <- combn(groups, 2, simplify = FALSE)

results_list <- list()
for(i in seq_along(comparisons)) {
  result <- hutcheson_test(alpha_diversity, "Group", 
                           comparisons[[i]][1], comparisons[[i]][2])
  results_list[[i]] <- result
}

results_df <- do.call(rbind, results_list)
print(results_df)

# 3. Compare specific combinations (e.g., BF Time1 vs BF Time2)
cat("\n--- Within-Group Temporal Comparisons ---\n")

for(grp in groups) {
  cat("\n", grp, ":\n")
  subset_data <- alpha_diversity[alpha_diversity$Group == grp, ]
  if(length(unique(subset_data$Time)) == 2) {
    result <- hutcheson_test(subset_data, "Time", "Time1", "Time2")
    print(result)
  }
}

# ============================================================
# Method 2: Using a more precise Hutcheson implementation
# ============================================================

# More accurate Hutcheson test that accounts for within-sample variance
hutcheson_precise <- function(abundance_matrix, group1_indices, group2_indices) {
  
  # Calculate Shannon and its variance for each sample
  calc_shannon_variance <- function(abundances) {
    # Remove zeros
    abundances <- abundances[abundances > 0]
    N <- sum(abundances)
    pi <- abundances / N
    
    # Shannon index
    H <- -sum(pi * log(pi))
    
    # Variance of Shannon (Hutcheson 1970 formula)
    var_H <- (sum(pi * (log(pi))^2) - H^2) / N + (length(pi) - 1) / (2 * N^2)
    
    return(c(H = H, var_H = var_H))
  }
  
  # Calculate for all group 1 samples
  group1_stats <- t(apply(abundance_matrix[group1_indices, ], 1, calc_shannon_variance))
  H1 <- mean(group1_stats[, "H"])
  var_H1 <- mean(group1_stats[, "var_H"])
  
  # Calculate for all group 2 samples
  group2_stats <- t(apply(abundance_matrix[group2_indices, ], 1, calc_shannon_variance))
  H2 <- mean(group2_stats[, "H"])
  var_H2 <- mean(group2_stats[, "var_H"])
  
  # T-statistic
  t_stat <- abs(H1 - H2) / sqrt(var_H1 + var_H2)
  
  # Degrees of freedom
  df <- (var_H1 + var_H2)^2 / 
    ((var_H1^2 / length(group1_indices)) + (var_H2^2 / length(group2_indices)))
  
  # P-value
  p_value <- 2 * (1 - pt(t_stat, df))
  
  return(list(
    H1 = H1,
    H2 = H2,
    t_statistic = t_stat,
    df = df,
    p_value = p_value
  ))
}

# Example usage with precise method
cat("\n=== Precise Hutcheson Test (Time 1 vs Time 2) ===\n")
time1_idx <- which(metadata$Time == "Time1")
time2_idx <- which(metadata$Time == "Time2")

precise_result <- hutcheson_precise(abundance_data, time1_idx, time2_idx)
cat("Mean Shannon Time1:", round(precise_result$H1, 4), "\n")
cat("Mean Shannon Time2:", round(precise_result$H2, 4), "\n")
cat("t-statistic:", round(precise_result$t_statistic, 4), "\n")
cat("df:", round(precise_result$df, 2), "\n")
cat("p-value:", round(precise_result$p_value, 4), "\n")
cat("Significant:", ifelse(precise_result$p_value < 0.05, "Yes", "No"), "\n")

# ============================================================
# Visualization of Results
# ============================================================

# Add significance annotations to your Shannon plot
library(ggpubr)

p_shannon_annotated <- ggplot(alpha_diversity, aes(x = Group, y = Shannon, fill = Time)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, 
               position = position_dodge(width = 0.75)) +
  geom_point(aes(shape = Time), 
             position = position_jitterdodge(
               jitter.width = 0.05,
               jitter.height = 0,
               dodge.width = 0.75
             ),
             size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("Time1" = "#66C2A5", "Time2" = "#FC8D62"),
                    labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_shape_manual(values = c("Time1" = 16, "Time2" = 17),
                     labels = c("Time1" = "Time 1", "Time2" = "Time 2")) +
  scale_x_discrete(labels = c("B" = "BF", "K" = "KF", "RM" = "RM", "WM" = "WM")) +
  theme_minimal() +
  labs(
    title = "Shannon Diversity Index",
    x = "Sample Type",
    y = "Shannon Index",
    fill = "Time Point",
    shape = "Time Point"
  ) +
  theme(legend.position = "right")

# Add significance stars if p < 0.05
if(hutcheson_time$p_value < 0.05) {
  p_shannon_annotated <- p_shannon_annotated +
    annotate("text", x = 2.5, y = max(alpha_diversity$Shannon) * 1.1,
             label = paste0("Time comparison: p = ", 
                            round(hutcheson_time$p_value, 4), " *"),
             size = 4, fontface = "italic")
}

print(p_shannon_annotated)
