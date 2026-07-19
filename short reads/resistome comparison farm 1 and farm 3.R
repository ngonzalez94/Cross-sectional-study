library(dplyr)
library(tidyr)
library(openxlsx)
library(ggVennDiagram)   
library(ggplot2)
library(patchwork)       
# install.packages(c("dplyr","tidyr","openxlsx","ggVennDiagram","ggplot2","patchwork"))

# ── 1. Load data ──────────────────────────────────────────────────────────────
data <- read.csv("resistome.csv", check.names = FALSE)
# If Excel:
# library(readxl); data <- read_excel("resistome.xlsx")

# Expected columns: Sample, Gene_Name, Read_Counts, Drug_Class, Category

# ── 2. Farm / sample definitions ─────────────────────────────────────────────
farms <- list(
  "Farm1_T1" = list(cmr = "E01", wm = "E02"),
  "Farm3_T1" = list(cmr = "E09", wm = "E10"),
  "Farm1_T2" = list(cmr = "F11", wm = "F12"),
  "Farm3_T2" = list(cmr = "G08", wm = "G09")
)

# ── 3. Helpers ────────────────────────────────────────────────────────────────
get_counts <- function(data, sample_id) {
  data %>%
    filter(Sample == sample_id) %>%
    select(Gene_Name, Read_Counts, Drug_Class) %>%
    group_by(Gene_Name, Drug_Class) %>%
    summarise(Read_Counts = sum(Read_Counts), .groups = "drop")
}

# ── 4. Build abundance tables per farm ───────────────────────────────────────
all_abund  <- list()   # for Excel
venn_sets  <- list()   # for Venn diagrams
heatmap_df <- list()   # for heatmap

for (farm in names(farms)) {
  cmr_id <- farms[[farm]]$cmr
  wm_id  <- farms[[farm]]$wm
  
  cmr_counts <- get_counts(data, cmr_id) %>% rename(CMR_Fed = Read_Counts)
  wm_counts  <- get_counts(data, wm_id)  %>% rename(WM_Fed  = Read_Counts)
  
  # Full join so all genes from both samples are included
  joined <- full_join(cmr_counts, wm_counts, by = c("Gene_Name", "Drug_Class")) %>%
    mutate(
      CMR_Fed  = replace_na(CMR_Fed, 0),
      WM_Fed   = replace_na(WM_Fed,  0),
      Diff_WM_minus_CMR = WM_Fed - CMR_Fed,
      Status = case_when(
        CMR_Fed == 0 ~ "WM Fed only",
        WM_Fed  == 0 ~ "CMR Fed only",
        TRUE         ~ "Shared"
      )
    ) %>%
    arrange(desc(abs(Diff_WM_minus_CMR)))
  
  all_abund[[farm]] <- joined
  
  # Venn sets (presence/absence)
  venn_sets[[farm]] <- list(
    "WM Fed"  = unique(wm_counts$Gene_Name),
    "CMR Fed" = unique(cmr_counts$Gene_Name)
  )
  
  # Heatmap data: keep diff column labelled by farm
  heatmap_df[[farm]] <- joined %>%
    select(Gene_Name, Diff_WM_minus_CMR) %>%
    rename(!!farm := Diff_WM_minus_CMR)
  
  cat(sprintf("\n%s | WM-only: %d | CMR-only: %d | Shared: %d\n",
              farm,
              sum(joined$Status == "WM Fed only"),
              sum(joined$Status == "CMR Fed only"),
              sum(joined$Status == "Shared")))
}

# ── 5. Excel output ───────────────────────────────────────────────────────────
wb         <- createWorkbook()
headerSty  <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
greenSty   <- createStyle(fgFill = "#C6EFCE")
redSty     <- createStyle(fgFill = "#FFC7CE")
blueSty    <- createStyle(fgFill = "#DDEBF7")

# Summary sheet
summary_df <- data.frame(
  Farm        = names(farms),
  CMR_Sample  = sapply(names(farms), function(f) farms[[f]]$cmr),
  WM_Sample   = sapply(names(farms), function(f) farms[[f]]$wm),
  WM_Only     = sapply(names(farms), function(f) sum(all_abund[[f]]$Status == "WM Fed only")),
  CMR_Only    = sapply(names(farms), function(f) sum(all_abund[[f]]$Status == "CMR Fed only")),
  Shared      = sapply(names(farms), function(f) sum(all_abund[[f]]$Status == "Shared")),
  Total_Genes = sapply(names(farms), function(f) nrow(all_abund[[f]]))
)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
addStyle(wb, "Summary", headerSty, rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
setColWidths(wb, "Summary", cols = 1:ncol(summary_df), widths = "auto")

# One sheet per farm
for (farm in names(farms)) {
  sname <- gsub("_", " ", farm)
  df    <- all_abund[[farm]]
  
  addWorksheet(wb, sname)
  writeData(wb, sname, df)
  addStyle(wb, sname, headerSty, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  setColWidths(wb, sname, cols = 1:ncol(df), widths = "auto")
  
  # Colour rows by status
  for (i in seq_len(nrow(df))) {
    sty <- if (df$Status[i] == "WM Fed only")  greenSty else
      if (df$Status[i] == "CMR Fed only") redSty   else blueSty
    addStyle(wb, sname, sty, rows = i + 1, cols = 1:ncol(df),
             gridExpand = TRUE, stack = TRUE)
  }
  
  # Conditional formatting on diff column: positive = green, negative = red
  diff_col <- which(names(df) == "Diff_WM_minus_CMR")
  conditionalFormatting(wb, sname,
                        cols = diff_col, rows = 2:(nrow(df) + 1),
                        rule = ">0", style = createStyle(fontColour = "#006100"))
  conditionalFormatting(wb, sname,
                        cols = diff_col, rows = 2:(nrow(df) + 1),
                        rule = "<0", style = createStyle(fontColour = "#9C0006"))
}

saveWorkbook(wb, "resistome_abundance.xlsx", overwrite = TRUE)
cat("\nExcel saved: resistome_abundance.xlsx\n")

# ── 6. Bar charts: top 15 genes by |diff| per farm (2x2 panel) ───────────────
bar_plots <- list()

for (farm in names(farms)) {
  df <- all_abund[[farm]] %>%
    mutate(Gene_short = sub(".*\\|", "", Gene_Name)) %>%
    group_by(Gene_short) %>%
    summarise(Diff_WM_minus_CMR = sum(Diff_WM_minus_CMR), .groups = "drop") %>%
    slice_max(abs(Diff_WM_minus_CMR), n = 15, with_ties = FALSE) %>%
    mutate(Direction = ifelse(Diff_WM_minus_CMR >= 0, "WM Fed higher", "CMR Fed higher"))
  
  p <- ggplot(df, aes(x = reorder(Gene_short, Diff_WM_minus_CMR),
                      y = Diff_WM_minus_CMR,
                      fill = Direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("WM Fed higher" = "#2ca25f", "CMR Fed higher" = "#de2d26")) +
    coord_flip() +
    labs(
      title = gsub("_", " ", farm),
      x = NULL,
      y = "Read Counts (WM Fed − CMR Fed)",
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  bar_plots[[farm]] <- p
}

bar_panel <- (bar_plots[[1]] | bar_plots[[2]]) /
  (bar_plots[[3]] | bar_plots[[4]]) +
  plot_annotation(
    title    = "Top 15 Genes by Abundance Difference (WM Fed − CMR Fed)",
    subtitle = "Green = WM Fed higher  |  Red = CMR Fed higher",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave("resistome_barcharts.png", bar_panel,
       width = 16, height = 12, dpi = 300, bg = "white")
cat("Bar charts saved: resistome_barcharts.png\n")

# ── 7. Heatmap: all genes across all 4 farms ──────────────────────────────────
heat_wide <- reduce(heatmap_df, full_join, by = "Gene_Name")
heat_wide <- heat_wide %>%
  mutate(across(-Gene_Name, ~replace_na(., 0))) %>%
  mutate(Gene_short = sub(".*\\|", "", Gene_Name)) %>%
  filter(rowSums(across(all_of(names(farms)), abs)) > 0)

# Pivot to long for ggplot
heat_long <- heat_wide %>%
  pivot_longer(cols = all_of(names(farms)),
               names_to = "Farm", values_to = "Diff") %>%
  mutate(Farm = gsub("_", " ", Farm))

# Optionally limit to top N genes by max absolute diff across farms
top_genes <- heat_wide %>%
  mutate(max_abs = apply(across(all_of(names(farms))), 1, function(x) max(abs(x)))) %>%
  slice_max(max_abs, n = 40) %>%
  pull(Gene_short)

heat_long_top <- heat_long %>% filter(Gene_short %in% top_genes)

heatmap_plot <- ggplot(heat_long_top,
                       aes(x = Farm, y = reorder(Gene_short, Diff), fill = Diff)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#de2d26",
    mid      = "white",
    high     = "#2ca25f",
    midpoint = 0,
    name     = "WM − CMR\nRead Counts"
  ) +
  labs(
    title    = "Resistome Abundance Difference Across Farms",
    subtitle = "Green = WM Fed higher  |  Red = CMR Fed higher  |  Top 40 genes by max |diff|",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey40"),
    axis.text.x     = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y     = element_text(size = 8),
    legend.position = "right"
  )

ggsave("resistome_heatmap.png", heatmap_plot,
       width = 10, height = 14, dpi = 300, bg = "white")
cat("Heatmap saved: resistome_heatmap.png\n")

# ── 8. Venn diagrams (2x2 panel) ─────────────────────────────────────────────
venn_plots <- list()

for (farm in names(farms)) {
  p <- ggVennDiagram(venn_sets[[farm]], label_alpha = 0, edge_size = 0.8) +
    scale_fill_gradient(low = "#f7fbff", high = "#2171b5") +
    scale_color_manual(values = c("WM Fed" = "#2ca25f", "CMR Fed" = "#de2d26")) +
    labs(title = gsub("_", " ", farm)) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "none"
    )
  venn_plots[[farm]] <- p
}

venn_panel <- (venn_plots[[1]] | venn_plots[[2]]) /
  (venn_plots[[3]] | venn_plots[[4]]) +
  plot_annotation(
    title    = "Resistome Gene Overlap: WM Fed vs CMR Fed",
    subtitle = "Numbers show unique and shared AMR genes per farm-timepoint",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave("resistome_venn_diagrams.png", venn_panel,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("Venn diagrams saved: resistome_venn_diagrams.png\n")