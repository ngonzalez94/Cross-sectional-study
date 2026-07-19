# Needed libraries
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(viridis)

# --- Load gene class info ---
gene_class_map <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/ARG_contigs/ARG_list_mech.xlsx", sheet = 1)

# --- Load display name map (caps -> real names) ---
gene_name_map <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/ARG_contigs/ARG_display_names.xlsx")

# --- File processor --- (only keeps rows where PRODUCT starts with "Drugs:")
process_files <- function(file_paths, source_label) {
  map_df(file_paths, ~ {
    read_csv(.x) %>%
      filter(str_starts(PRODUCT, "Drugs:")) %>%
      mutate(Source = source_label, File = basename(.x))
  })
}

# --- Load data ---
keeper_df <- process_files(list.files(
  "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/ARG_contigs/Merged_output/Faeces/Keeper_calves/MEGARES", 
  pattern = "_ARGs_with_families\\.csv$", full.names = TRUE), 
  "Replacement calves' faeces")

bobby_df <- process_files(list.files(
  "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/ARG_contigs/Merged_output/Faeces/Bobby_calves/MEGARES", 
  pattern = "_ARGs_with_families\\.csv$", full.names = TRUE), 
  "Non-replacement calves' faeces")

rm_df <- process_files(list.files(
  "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/ARG_contigs/Merged_output/Milks/RM/MEGARES", 
  pattern = "_ARGs_with_families\\.csv$", full.names = TRUE), 
  "Replacement Milk")

wm_df <- process_files(list.files(
  "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/ARG_contigs/Merged_output/Milks/WM/MEGARES", 
  pattern = "_ARGs_with_families\\.csv$", full.names = TRUE), 
  "Waste Milk")

# --- Manually define which families to show per source ---
top10_per_source <- tribble(
  ~Source,             ~Family,
  "Replacement calves' faeces",     "Enterobacteriaceae",
  "Replacement calves' faeces",     "Bacteroidaceae",
  "Replacement calves' faeces",     "Oscillospiraceae",
  "Replacement calves' faeces",     "Staphylococcaceae",
  "Replacement calves' faeces",     "Streptococcaceae",
  "Replacement calves' faeces",     "Lachnospiraceae",
  "Replacement calves' faeces",     "Fibrobacteraceae",
  "Replacement calves' faeces",     "Candidatus Cloacimonadaceae", 
  "Replacement calves' faeces",     "Gemmatimonadaceae",
  "Replacement calves' faeces",     "Clostridiaceae", 
  "Non-replacement calves' faeces",      "Enterobacteriaceae",
  "Non-replacement calves' faeces",      "Bacteroidaceae",
  "Non-replacement calves' faeces",      "Oscillospiraceae",
  "Non-replacement calves' faeces",      "Streptococcaceae",
  "Non-replacement calves' faeces",      "Staphylococcaceae",
  "Non-replacement calves' faeces",      "Enterococcaceae",
  "Non-replacement calves' faeces",      "Lachnospiraceae",
  "Non-replacement calves' faeces",      "Candidatus Cloacimonadaceae",
  "Non-replacement calves' faeces",      "Clostridiaceae", 
  "Non-replacement calves' faeces",      "Gemmatimonadaceae",
  "Replacement Milk",  "Enterobacteriaceae",
  "Replacement Milk",  "Streptococcaceae",
  "Replacement Milk",  "Enterococcaceae",
  "Replacement Milk",  "Hafniaceae",
  "Replacement Milk",  "Staphylococcaceae",
  "Replacement Milk",  "Moraxellaceae",
  "Replacement Milk",  "Clostridiaceae",
  "Replacement Milk",  "Acidithiobacillaceae",
  "Replacement Milk",  "Rhizobiaceae",
  "Waste Milk",        "Enterobacteriaceae",
  "Waste Milk",        "Streptococcaceae",
  "Waste Milk",        "Hafniaceae",
  "Waste Milk",        "Moraxellaceae",
  "Waste Milk",        "Lactobacillaceae",
  "Waste Milk",        "Chitinibacteraceae",
)

# --- All families appearing in any source's manual selection ---
all_selected_families <- distinct(top10_per_source, Family) %>% pull(Family)

# --- Combine, clean, filter — all in one block ---
all_data <- bind_rows(keeper_df, bobby_df, rm_df, wm_df) %>%
  mutate(Family = str_trim(Family)) %>%                        # fix whitespace
  filter(!str_ends(str_to_lower(Family), "viridae")) %>%       # remove viruses
  left_join(gene_class_map, by = "GENE") %>%                   # add class info
  filter(Family %in% all_selected_families)                    # keep only selected families


#Too look for genes and families
result <-all_data %>% filter(Family == "Enterobacteriaceae", GENE == "TETA")
write.csv(result, "filtered_contigs_Enterobacteriaceae_TETA.csv", row.names = FALSE)
# --- Factor formatting ---
all_data$Source <- factor(all_data$Source, 
                          levels = c("Replacement calves' faeces", "Non-replacement calves' faeces", "Replacement Milk", "Waste Milk"))

# --- Count gene occurrences (integer) ---
gene_counts <- all_data %>%
  group_by(GENE, Source, Family, Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Count = as.integer(Count))

# --- Apply display names from the name map ---
gene_counts <- gene_counts %>%
  left_join(gene_name_map, by = c("GENE" = "Name_og")) %>%
  mutate(Gene_label = if_else(is.na(Name_change), GENE, Name_change)) %>%
  select(-Name_change)

# --- Remove BIOCIDES and NA class ---
gene_counts <- gene_counts %>%
  filter(!Class %in% c("BIOCIDES"), !is.na(Class))%>%
  filter(!GENE %in% c("LMRD"), !is.na(GENE))

# --- Keep only rows where Family is in that source's manual selection ---
gene_counts <- gene_counts %>%
  semi_join(top10_per_source, by = c("Source", "Family"))

# --- Gene order: sort by Class then Gene_label ---
gene_counts <- gene_counts %>%
  arrange(Class, Gene_label) %>%
  mutate(Gene_label = factor(Gene_label, levels = unique(Gene_label)))

# --- Family order: ranked by total contigs, highest first (left to right) ---
#family_order <- gene_counts %>%
 # group_by(Family) %>%
#  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
#  arrange(desc(Total)) %>%
#  pull(Family)
family_order <- top10_per_source %>%
  distinct(Family) %>%
  pull(Family)

gene_counts$Family <- factor(gene_counts$Family, levels = family_order)

# --- Source factor ---
gene_counts$Source <- factor(gene_counts$Source, 
                             levels = c("Non-replacement calves' faeces", "Replacement calves' faeces", "Replacement Milk", "Waste Milk"))

# --- Plot ---
ggplot(gene_counts, aes(x = Family, y = Gene_label, fill = Count)) +
  geom_tile() +
  scale_fill_viridis(
    option    = "cividis",
    direction = -1,
    na.value  = "white",
    name      = "Contigs",
    breaks    = scales::breaks_pretty(),
    labels    = scales::label_number(accuracy = 1)
  ) +
  facet_grid(
    rows   = vars(Class),
    cols   = vars(Source),
    scales = "free",
    space  = "free",
    switch = "y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement.y  = "outside",
    strip.background   = element_blank(),
    strip.text.y.left  = element_text(face = "bold"),
    strip.text.x       = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.text.y        = element_text(face = "italic"),
    panel.grid         = element_blank(),
    axis.title.x       = element_text(margin = margin(t = 10)),
    axis.title.y       = element_text(margin = margin(r = 10))
  ) +
  labs(
    title = "ARG Gene Occurrence by Family, Class and Source",
    x     = "Taxonomic Family (ordered by total contigs)",
    y     = "ARG Gene (grouped by Class)"
  )

# --- Export summary table ---
summary_counts <- gene_counts %>%
  group_by(Source, Gene_label, GENE, Family, Class) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop")

print(summary_counts)
write.csv(summary_counts, "summary_counts_arg_contig.csv", row.names = FALSE)

# --- Count total contigs per Family and Source ---
family_source_counts <- gene_counts %>%
  group_by(Family, Source) %>%
  summarise(Total_Contigs = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  arrange(Source, desc(Total_Contigs))

print(family_source_counts)
write.csv(family_source_counts, "family_source_contig_counts.csv", row.names = FALSE)

# --- Total contigs with a family classification (excluding BIOCIDES and NA) ---
total_classified <- all_data %>%
  filter(!is.na(Family), !Class %in% c("BIOCIDES"), !is.na(Class)) %>%
  nrow()

cat("Total contigs with a family classification (excluding BIOCIDES and NA):", total_classified, "\n")
## -----------
family_source_counts <- gene_counts %>%
  group_by(Source) %>%
  summarise(Total_Contigs = sum(Count), .groups = "drop") %>%
  arrange(Source, desc(Total_Contigs))

print(family_source_counts)

# --- Total contigs with a family classification (excluding BIOCIDES and NA) ---
total_classified <- all_data %>%
  filter(!is.na(Family), !Class %in% c("BIOCIDES"), !is.na(Class)) %>%
  nrow()

cat("Total contigs with a family classification (excluding BIOCIDES and NA):", total_classified, "\n")

# --- Total contigs per source across ALL families (before any filtering) ---
all_source_summary <- bind_rows(keeper_df, bobby_df, rm_df, wm_df) %>%
  mutate(Family = str_trim(Family)) %>%
  filter(!str_ends(str_to_lower(Family), "viridae")) %>%
  filter(str_starts(PRODUCT, "Drugs:")) %>%
  filter(!is.na(Family)) %>%
  group_by(Source) %>%
  summarise(Total_Contigs = n(), .groups = "drop") %>%
  arrange(factor(Source, levels = c("Bobby faeces", "Keeper faeces", "Replacement Milk", "Waste Milk")))

print(all_source_summary)
write.csv(all_source_summary, "all_source_contig_summary.csv", row.names = FALSE)

# --- Total unique families per source (all families, no filtering) ---
all_family_summary <- bind_rows(keeper_df, bobby_df, rm_df, wm_df) %>%
  mutate(Family = str_trim(Family)) %>%
  filter(!str_ends(str_to_lower(Family), "viridae")) %>%
  filter(str_starts(PRODUCT, "Drugs:")) %>%
  filter(!is.na(Family)) %>%
  group_by(Source) %>%
  summarise(Total_Families = n_distinct(Family), .groups = "drop") %>%
  arrange(factor(Source, levels = c("Bobby faeces", "Keeper faeces", "Replacement Milk", "Waste Milk")))

print(all_family_summary)
write.csv(all_family_summary, "all_family_summary.csv", row.names = FALSE)

# --- Total unique families across all sources combined ---
total_unique_families <- bind_rows(keeper_df, bobby_df, rm_df, wm_df) %>%
  mutate(Family = str_trim(Family)) %>%
  filter(!str_ends(str_to_lower(Family), "viridae")) %>%
  filter(str_starts(PRODUCT, "Drugs:")) %>%
  filter(!is.na(Family)) %>%
  distinct(Family) %>%
  nrow()

cat("Total unique families across all sources:", total_unique_families, "\n")
