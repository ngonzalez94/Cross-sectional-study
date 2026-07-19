library(tidyverse)

setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/Taxa_L_S")

files <- list.files(pattern = "_family_reads.tsv$", full.names = TRUE)

df <- map_df(files, function(f) {
  read_tsv(
    f,
    col_names = c("Family", "Reads"),
    col_types = cols()
  ) %>%
    mutate(File = basename(f))
})
colnames(df)
head(df)

df <- df %>%
  mutate(
    Sample = case_when(
      str_detect(File, "MD0117") ~ "Farm2_T1_LR",
      str_detect(File, "MD0145") ~ "Farm3_T2_LR",
      str_detect(File, "E06")    ~ "Farm2_T1_SR",
      str_detect(File, "G09")    ~ "Farm3_T2_SR",
      TRUE ~ File
    ),
    Technology = case_when(
      str_detect(File, "^MD")   ~ "LR",
      str_detect(File, "^trim") ~ "SR",
      TRUE ~ "Other"
    )
  )


df_rel <- df %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Reads / sum(Reads)) %>%
  ungroup()


top10 <- df_rel %>%
  group_by(Family) %>%
  summarise(TotalRA = sum(RelAbundance)) %>%
  arrange(desc(TotalRA)) %>%
  slice_head(n = 10) %>%
  pull(Family)

df_plot <- df_rel %>%
  mutate(Family_plot = if_else(Family %in% top10, Family, "Other")) %>%
  group_by(Sample, Technology, Family_plot) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups = "drop")

# Named vector of colors
top_colors <- RColorBrewer::brewer.pal(10, "Set3")
names(top_colors) <- setdiff(top10, "Other")  # match the factor levels
top_colors <- c(top_colors, Other = "grey70") # grey for Other


# Keep top10 in order, put "Other" last
df_plot$Family_plot <- factor(
  df_plot$Family_plot,
  levels = c(setdiff(top10, "Other"), "Other")
)


df_plot$Sample <- factor(
  df_plot$Sample,
  levels = c(
    "Farm2_T1_LR",
    "Farm2_T1_SR",
    "Farm3_T2_LR",
    "Farm3_T2_SR"
  )
)

df_plot <- df_plot %>%
  mutate(RelAbundance = RelAbundance * 100)


ggplot(df_plot, aes(x = Sample, y = RelAbundance, fill = Family_plot)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = top_colors) +
  labs(x = "", y = "Relative abundance", fill = "Family") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),    # Legend title
    legend.text = element_text(size = 14),
    guide_legend(label.theme = element_text(face = "italic")),
    panel.grid = element_blank()
  )

library(scales)

ggplot(df_plot, aes(x = Sample, y = RelAbundance, fill = Family_plot)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = top_colors) +
  scale_y_continuous(labels = function(x) paste0(round(x,1), "%")) +  # add % to y-axis
  labs(x = "", y = "Relative abundance (%)", fill = "Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




library(tidyverse)
library(tidyr)

# If you want relative abundance
df_wide <- df_plot %>%
  select(Sample, Family_plot, RelAbundance) %>%
  pivot_wider(
    names_from = Family_plot,
    values_from = RelAbundance,
    values_fill = 0  # fills missing families with 0
  )

# Pivot so families = columns, samples = rows
df_table <- df_plot %>%
  select(Sample, Family_plot, RelAbundance) %>%
  pivot_wider(
    names_from = Sample,     # columns = samples
    values_from = RelAbundance,
    values_fill = 0
  )

df_table
