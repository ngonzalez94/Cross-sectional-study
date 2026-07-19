library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── 1. Load data ──────────────────────────────────────────────────────────────
# Replace with your actual file paths
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken")
keepers <- read.csv("keepers_15Families.csv", check.names = FALSE)
bobbies <- read.csv("bobby_15Families.csv",   check.names = FALSE)

# Drop row-index column if present
if (names(keepers)[1] %in% c("", "X")) keepers <- keepers[, -1]
if (names(bobbies)[1] %in% c("", "X")) bobbies <- bobbies[, -1]

# ── 2. Comparisons of interest ────────────────────────────────────────────────
comparisons <- c("Farm1_T1", "Farm3_T1", "Farm1_T2", "Farm3_T2")

# ── 3. Find common taxa ───────────────────────────────────────────────────────
common_taxa <- intersect(keepers$Taxon, bobbies$Taxon)

keepers_sub <- keepers %>% filter(Taxon %in% common_taxa)
bobbies_sub <- bobbies %>% filter(Taxon %in% common_taxa)

# ── 4. Build difference table ─────────────────────────────────────────────────
# For each farm column: WM Fed (keeper) minus CMR Fed (bobby)
diff_df <- data.frame(Taxon = common_taxa)

for (col in comparisons) {
  k_vals <- keepers_sub[[col]][match(common_taxa, keepers_sub$Taxon)]
  b_vals <- bobbies_sub[[col]][match(common_taxa, bobbies_sub$Taxon)]
  diff_df[[col]] <- k_vals - b_vals
}

# ── 5. Bar charts: top 10 families per farm (2x2 panel) ───────────────────────
bar_plots <- list()

for (col in comparisons) {
  df <- diff_df %>%
    select(Taxon, Diff = all_of(col)) %>%
    slice_max(abs(Diff), n = 15, with_ties = FALSE) %>%
    mutate(Direction = ifelse(Diff >= 0, "WM Fed higher", "CMR Fed higher"))
  
  p <- ggplot(df, aes(x = reorder(Taxon, Diff),
                      y = Diff,
                      fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = c("WM Fed higher"  = "#2ca25f",
                                 "CMR Fed higher" = "#de2d26"),
                      drop = FALSE) +
    coord_flip() +
    labs(
      title = gsub("_", " ", col),
      x     = NULL,
      y     = "Relative Abundance (WM Fed − CMR Fed)",
      fill  = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      legend.position    = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  bar_plots[[col]] <- p
}

# ── 6. Combine into 2x2 panel ─────────────────────────────────────────────────
combined <- (bar_plots[[1]] | bar_plots[[2]]) /
  (bar_plots[[3]] | bar_plots[[4]]) +
  plot_annotation(
    title    = "Top 10 Families by Abundance Difference (WM Fed − CMR Fed)",
    subtitle = "Green = WM Fed higher  |  Red = CMR Fed higher",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave("taxa_barcharts.png", combined,
       width = 14, height = 10, dpi = 300, bg = "white")
cat("Saved: taxa_barcharts.png\n")