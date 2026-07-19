# Load required libraries
library(tidyverse)
library(readxl)
library(ggplot2)
install.packages("gghalves")
library(gghalves)
# Load necessary libraries
library(ggplot2)
install.packages("gghalves", repos = c('https://erocoar.r-universe.dev', 'https://cloud.r-project.org'))

# Define the file path
file_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/bracken_all_output.xlsx"
# Load the bracken data
# No 'header' or 'stringsAsFactors' arguments
data <- read_excel(file_path, sheet=8)
View(data)

# Set the first column as row names (family names) and keep only numeric values
# Store taxon names first
taxon_names <- data[[1]]

# Remove the first column and set taxon names as row names
data <- data[,-1]
rownames(data) <- taxon_names


# Transpose the data so that families are columns and samples are rows
data <- t(data)

# Convert to a data frame and restore taxon names
data <- as.data.frame(data)
data$Sample <- rownames(data)  # Restore sample names

# Reshape data to long format
data_long <- data %>%
  pivot_longer(
    cols = -Sample,  # Keep sample column, pivot everything else
    names_to = "Taxon",
    values_to = "Abundance"
  )

#Group your samples
group_column <- c("Bobby", "Keeper", "Bobby", "Keeper", "Bobby", "Keeper", "Bobby", "Keeper", "Bobby" ,"Keeper","Bobby" ,"Keeper", "Bobby", "Keeper", "Bobby", "Keeper", "Bobby", "Keeper", "Bobby", "Keeper", "Bobby" ,"Keeper" )  
# Add the Group column at the beginning of the data
data <- cbind(Group = group_column, data)
# Convert to a data frame for better handling
data <- as.data.frame(data)

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(
    cols = -c(Group,Sample), # Exclude 'Group' column from reshaping
    names_to = "Taxon", # New column for taxon names
    values_to = "Abundance" # New column for abundance values
  )

# Ensure the Abundance column is numeric
data_long <- data_long %>%
  mutate(Abundance = as.numeric(Abundance))

# Calculate the total abundance for each taxon across all samples
top_taxa <- data_long %>%
  group_by(Taxon) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%  # Keep only the top 10 taxa
  pull(Taxon)  # Extract the taxon names

# Filter data to keep only the top 10 taxa and group the rest as "Other"
data_long <- data_long %>%
  mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
  group_by(Sample, Group, Taxon) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")  # Aggregate the "Other" group

# Order taxa: "Other" at the top, remaining taxa by total abundance
taxa_order <- c("Other", top_taxa)  # "Other" is first
data_long$Taxon <- factor(data_long$Taxon, levels = taxa_order)

# Plot the stacked bar chart
ggplot(data_long, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_fill_manual(values = custom_colors) +
  labs(
    title = "Stacked Taxa Bar Chart (Top 10 Taxa) Faecal samples",
    x = "Sample",
    y = "Total Reads",
    fill = "Taxon"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("grey80",'#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd', '#c51b7d',"#f1b6da"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X-axis labels
  facet_wrap(~ Group, scales = "free_x", ncol = 2)  # Split into facets by 'Group'


#Parse the sample name for the taxonomical distribution
df_tidy <- data_long %>%
  separate(
    Sample,
    into = c("Farm", "Group", "Time"),
    sep = "_"
  )
#Plot
ggplot(df_tidy, aes(x = Group, y = Abundance, fill = Group)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.6
  )+
  facet_grid(
    Taxon ~ Time,
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c("B" = "darkorange", "K" = "darkred")
  ) +
  labs(
    x = "",
    y = "Taxonomic abundance"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

#Now lets add paired lines

ggplot(df_tidy, aes(x = Group, y = Abundance)) +
  
  # Violins
  geom_violin(
    aes(fill = Group),
    trim = FALSE,
    alpha = 0.7
  ) +
  
  # Paired lines per farm
  geom_line(
    aes(group = Farm),
    color = "grey40",
    alpha = 0.6,
    linewidth = 0.4
  ) +
  
  facet_grid(
    Taxon ~ Time,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c("B" = "darkorange", "K" = "darkred")
  ) +
  
  labs(
    x = "",
    y = "Taxonomic abundance"
  ) +
  
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )

#Better try to change to a beeswarm and change the colours
#install.packages("ggbeeswarm")
library(ggbeeswarm)

#let's plot them
ggplot(df_tidy, aes(x = Group, y = Abundance)) +
  
  # Beeswarm points
  geom_quasirandom(
    aes(color = Farm),
    width = 0.25,
    size = 2,
    alpha = 0.85
  ) +
  
  facet_grid(
    Taxon ~ Time,
    scales = "free_y"
  ) +
  
  labs(
    x = "",
    y = "Taxonomic abundance",
    color = "Farm"
  ) +
  
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )

#simpler but useless
ggplot(df_tidy, aes(x = Group, y = Abundance)) +
  geom_quasirandom(
    aes(color = Farm),
    width = 0.2,
    size = 1.4,
    alpha = 0.7
  )
theme(legend.position = "none")


#Let's try and add all the plots together

ggplot(df_tidy, aes(x = Group, y = Abundance)) +
  
  # 1. Violin (distribution)
  geom_violin(
    aes(fill = Group),
    trim = FALSE,
    alpha = 0.35,
    width = 0.9
  ) +
  
  # 2. Boxplot (summary)
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  
  # 3. Paired lines per farm (optional but works well)
  #geom_line(
  #aes(group = Farm),
  #color = "grey50",
  #alpha = 0.5,
  # linewidth = 0.4
  #  ) +
  
  # 4. Beeswarm (raw data)
  geom_quasirandom(
    aes(color = Farm),
    width = 0.25,
    size = 2,
    alpha = 0.85
  ) +
  
  facet_grid(
    Taxon ~ Time,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c("B" = "orange", "K" = "red")
  ) +
  
  scale_color_viridis_d() +
  
  labs(
    x = "",
    y = "Taxonomic abundance",
    fill = "Group",
    color = "Farm"
  ) +
  
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )

##What if I put them next to eachother?

df_tidy <- df_tidy %>%
  mutate(
    Group_Time = interaction(Group, Time, sep = "_")
  )

df_tidy$Group_Time <- factor(
  df_tidy$Group_Time,
  levels = c("B_T1", "B_T2", "K_T1", "K_T2")
)

#Plot
ggplot(df_tidy, aes(x = Group_Time, y = Abundance)) +
  
  # Violin
  geom_violin(
    aes(fill = Group),
    trim = FALSE,
    alpha = 0.35,
    width = 0.9
  ) +
  
  # Boxplot
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  
  # Paired lines per farm
  #geom_line(
  # aes(group = Farm),
  # color = "grey50",
  #alpha = 0.5,
  #linewidth = 0.4
  # ) +
  
  # Beeswarm
  geom_quasirandom(
    aes(color = Farm),
    width = 0.25,
    size = 2,
    alpha = 0.85
  ) +
  
  facet_grid(
    Taxon ~ .,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c("B" = "orange", "K" = "red")
  ) +
  
  scale_color_viridis_d() +
  
  labs(
    x = "",
    y = "Taxonomic abundance"
  ) +
  
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
##Let's test for statistical differences 
library(dplyr)
library(tidyr)

median_tests <- df_tidy %>%
  select(Farm, Taxon, Group, Time, Abundance) %>%
  
  pivot_wider(
    names_from = Time,
    values_from = Abundance
  ) %>%
  
  group_by(Taxon, Group) %>%
  summarise(
    n_pairs   = sum(!is.na(T1) & !is.na(T2)),
    
    median_T1 = median(T1, na.rm = TRUE),
    IQR_T1    = IQR(T1, na.rm = TRUE),
    
    median_T2 = median(T2, na.rm = TRUE),
    IQR_T2    = IQR(T2, na.rm = TRUE),
    
    p_value   = wilcox.test(T1, T2, paired = TRUE)$p.value,
    
    .groups = "drop"
  ) %>%
  
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )
View(median_tests)

# Compare between farms at each time point
# Compare among farms within each Group (pooling across time)
among_farms_tests <- df_tidy %>%
  select(Farm, Taxon, Group, Abundance) %>%
  
  group_by(Taxon, Group) %>%
  summarise(
    n_farms = n_distinct(Farm),
    n_samples = n(),
    
    median_abundance = median(Abundance, na.rm = TRUE),
    IQR_abundance = IQR(Abundance, na.rm = TRUE),
    
    p_value = if(n_distinct(Farm) > 1) {
      kruskal.test(Abundance ~ Farm)$p.value
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

View(among_farms_tests)

#There's no diference, therefore we can move on and pool them together, then compare farm to farm.

df_BK <- df_tidy %>%
  filter(Group %in% c("B", "K")) %>%
  mutate(
    Group_factor = factor(Group, levels = c("B", "K")),
    # Assign shapes: T1 = triangle, T2 = circle
    point_shape = ifelse(Time == "T1", 24, 21)  # 24 = triangle up, 21 = circle
  )


#Plot
library(ggplot2)
library(ggbeeswarm)

df_BK$Group[df_BK$Group == "B"] <- "Non-replacement calves' faeces"
df_BK$Group[df_BK$Group == "K"] <- "Replacement calves' faeces"
df_BK$Group_factor <- as.character(df_BK$Group_factor)
df_BK$Group_factor[df_BK$Group_factor == "B"] <- "Non-RF"
df_BK$Group_factor[df_BK$Group_factor == "K"] <- "RF"



ggplot(df_BK, aes(x = Group_factor, y = Abundance, fill = Group)) +
  
  # Violin per group
  geom_violin(trim = FALSE, alpha = 0.3, width = 0.9) +
  
  # Boxplot overlay
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.6) +
  
  # Paired lines per farm
  #geom_line(aes(group = Farm), color = "grey50", alpha = 0.5, linewidth = 0.4) +
  
  # Beeswarm points: shape = T1/T2, color = farm
  geom_quasirandom(aes(color = Farm, shape = factor(point_shape)),
                   width = 0.2, size = 2, fill = "white") +
  
  # Map the shapes manually so the legend shows T1/T2 clearly
  scale_shape_manual(
    name = "Time Point",
    values = c(`24` = 24, `21` = 21),
    labels = c("T1", "T2")
  ) +
  
  # Colors for violins
  scale_fill_manual(values = c("Non-replacement calves' faeces" = "orange", "Replacement calves' faeces" = "red")) +
  scale_color_viridis_d() +
  
  # Facet per taxon
  facet_wrap(~ Taxon, scales = "free_y") +
  
  labs(x = "", y = " Relative abundance") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))

#test
# Keep only one abundance per farm × group × taxon (for this test)
df_BK_wide <- df_BK %>%
  group_by(Farm, Taxon, Group) %>%
  summarise(Abundance_total = sum(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Group,
    values_from = Abundance_total
  )

#Run paired Wilcoxon test per taxon
wilcox_BK <- df_BK_wide %>%
  group_by(Taxon) %>%
  summarise(
    n_pairs   = sum(!is.na(B) & !is.na(K)),
    median_B  = median(B, na.rm = TRUE),
    IQR_B     = IQR(B, na.rm = TRUE),
    median_K  = median(K, na.rm = TRUE),
    IQR_K     = IQR(K, na.rm = TRUE),
    p_value   = wilcox.test(B, K, paired = TRUE)$p.value,
    .groups   = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # optional BH correction

#Now let's check the effect size
# Prepare wide-format data: one row per farm × taxon, columns B and K
df_BK_wide <- df_BK %>%
  group_by(Farm, Taxon, Group) %>%
  summarise(Abundance_total = sum(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Group,
    values_from = Abundance_total
  )

# Run paired Wilcoxon and calculate effect size
wilcox_BK <- df_BK_wide %>%
  group_by(Taxon) %>%
  summarise(
    median_B  = median(B, na.rm = TRUE),
    IQR_B     = IQR(B, na.rm = TRUE),
    median_K  = median(K, na.rm = TRUE),
    IQR_K     = IQR(K, na.rm = TRUE),
    median_diff = median(B - K, na.rm = TRUE),        # Effect size
    p_value   = wilcox.test(B, K, paired = TRUE)$p.value,
    .groups   = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),        # optional multiple testing
    significant = ifelse(p_value < 0.05, "Yes", "No")
  )
## T1 and T2

# Filter for T1 only
df_T1 <- df_BK %>%
  filter(Time == "T1") %>%
  select(Farm, Group, Taxon, Abundance)

# Reshape to wide: one row per Farm × Taxon, columns Bobby and Keeper
df_T1_wide <- df_T1 %>%
  pivot_wider(
    names_from = Group,   # B vs K
    values_from = Abundance
  )
wilcox_T1 <- df_T1_wide %>%
  group_by(Taxon) %>%
  summarise(
    median_Bobby  = median(B, na.rm = TRUE),
    IQR_Bobby     = IQR(B, na.rm = TRUE),
    median_Keeper = median(K, na.rm = TRUE),
    IQR_Keeper    = IQR(K, na.rm = TRUE),
    median_diff   = median(B - K, na.rm = TRUE),
    p_value       = wilcox.test(B, K, paired = TRUE)$p.value,
    .groups       = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = ifelse(p_value < 0.05, "Yes", "No")
  )

#NOW T2
df_T2 <- df_BK %>%
  filter(Time == "T2") %>%
  select(Farm, Group, Taxon, Abundance)

df_T2_wide <- df_T2 %>%
  pivot_wider(names_from = Group, values_from = Abundance)

wilcox_T2 <- df_T2_wide %>%
  group_by(Taxon) %>%
  summarise(
    median_Bobby  = median(B, na.rm = TRUE),
    IQR_Bobby     = IQR(B, na.rm = TRUE),
    median_Keeper = median(K, na.rm = TRUE),
    IQR_Keeper    = IQR(K, na.rm = TRUE),
    median_diff   = median(B - K, na.rm = TRUE),
    p_value       = wilcox.test(B, K, paired = TRUE)$p.value,
    .groups       = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = ifelse(p_value < 0.05, "Yes", "No")
  )

#NOW T2
df_T2 <- df_BK %>%
  filter(Time == "T2") %>%
  select(Farm, Group, Taxon, Abundance)

df_T2_wide <- df_T2 %>%
  pivot_wider(names_from = Group, values_from = Abundance)

wilcox_T2 <- df_T2_wide %>%
  group_by(Taxon) %>%
  summarise(
    median_Bobby  = median(B, na.rm = TRUE),
    IQR_Bobby     = IQR(B, na.rm = TRUE),
    median_Keeper = median(K, na.rm = TRUE),
    IQR_Keeper    = IQR(K, na.rm = TRUE),
    median_diff   = median(B - K, na.rm = TRUE),
    p_value       = wilcox.test(B, K, paired = TRUE)$p.value,
    .groups       = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = ifelse(p_value < 0.05, "Yes", "No")
  )
