setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken/")

library(dplyr)
library(tidyr)
library(openxlsx)

# ── 1. Load data ──────────────────────────────────────────────────────────────
# Replace these paths with your actual file paths
keepers <- read.csv("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken/keeper_15Families.csv", check.names = FALSE)
bobbies <- read.csv("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Bracken/bobby_15Families.csv",   check.names = FALSE)

# Drop the row-index column if present (unnamed first column)
if (names(keepers)[1] == "" | names(keepers)[1] == "X") keepers <- keepers[, -1]
if (names(bobbies)[1] == "" | names(bobbies)[1] == "X") bobbies <- bobbies[, -1]

# ── 2. Define the four comparisons ───────────────────────────────────────────
comparisons <- c("Farm1_T1", "Farm3_T1", "Farm1_T2", "Farm3_T2")

# ── 3. Find taxa common to both datasets ─────────────────────────────────────
common_taxa <- intersect(keepers$Taxon, bobbies$Taxon)
cat("Common taxa:", length(common_taxa), "\n")
cat(common_taxa, sep = "\n")

keepers_sub <- keepers %>% filter(Taxon %in% common_taxa)
bobbies_sub <- bobbies %>% filter(Taxon %in% common_taxa)

# ── 4. Build difference tables ───────────────────────────────────────────────
results <- data.frame(Taxon = common_taxa)

for (col in comparisons) {
  
  k_vals <- keepers_sub[[col]][match(common_taxa, keepers_sub$Taxon)]
  b_vals <- bobbies_sub[[col]][match(common_taxa, bobbies_sub$Taxon)]
  
  results[[paste0(col, "_Keeper")]]   <- k_vals
  results[[paste0(col, "_Bobby")]]    <- b_vals
  results[[paste0(col, "_Diff_K-B")]] <- k_vals - b_vals          # signed
  results[[paste0(col, "_AbsDiff")]]  <- abs(k_vals - b_vals)     # absolute
}

# ── 5. Write to Excel ─────────────────────────────────────────────────────────
wb <- createWorkbook()

# --- Sheet 1: Full results table ---
addWorksheet(wb, "All_Comparisons")
writeData(wb, "All_Comparisons", results)

# Style: header row bold
headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
addStyle(wb, "All_Comparisons", headerStyle, rows = 1, cols = 1:ncol(results), gridExpand = TRUE)
setColWidths(wb, "All_Comparisons", cols = 1:ncol(results), widths = "auto")

# Colour-code signed diff columns: green = keeper higher, red = bobby higher
signedCols <- grep("_Diff_K-B", names(results))
for (j in signedCols) {
  for (i in seq_len(nrow(results))) {
    val <- results[i, j]
    fill <- if (!is.na(val) && val > 0) "#C6EFCE" else if (!is.na(val) && val < 0) "#FFC7CE" else "#FFFFFF"
    addStyle(wb, "All_Comparisons",
             createStyle(fgFill = fill),
             rows = i + 1, cols = j, stack = TRUE)
  }
}

# --- Sheets 2-5: One sheet per comparison ---
for (col in comparisons) {
  sheet_cols <- c("Taxon",
                  paste0(col, "_Keeper"),
                  paste0(col, "_Bobby"),
                  paste0(col, "_Diff_K-B"),
                  paste0(col, "_AbsDiff"))
  
  sheet_data <- results[, sheet_cols]
  sheet_name <- gsub("_", " ", col)   # e.g. "Farm1 T1"
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sheet_data)
  addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(sheet_data), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(sheet_data), widths = "auto")
  
  # Colour signed diff column on individual sheets
  diff_col_idx <- which(names(sheet_data) == paste0(col, "_Diff_K-B"))
  for (i in seq_len(nrow(sheet_data))) {
    val <- sheet_data[i, diff_col_idx]
    fill <- if (!is.na(val) && val > 0) "#C6EFCE" else if (!is.na(val) && val < 0) "#FFC7CE" else "#FFFFFF"
    addStyle(wb, sheet_name,
             createStyle(fgFill = fill),
             rows = i + 1, cols = diff_col_idx, stack = TRUE)
  }
}

# ── 6. Save ───────────────────────────────────────────────────────────────────
saveWorkbook(wb, "keepers_vs_bobbies_differences.xlsx", overwrite = TRUE)
cat("\nDone! File saved: keepers_vs_bobbies_differences.xlsx\n")