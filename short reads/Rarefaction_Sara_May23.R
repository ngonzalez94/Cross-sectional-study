#devtools::install_github("G-Thomson/Manu")
library("tidyverse")
library("ggplot2")
library("plyr")
library("readr")
library("shiny")
library("gghighlight")

library("dplyr")
install.packages("RColorBrewer")
library(RColorBrewer)

selected_colours <- brewer.pal(8, "Set2")[c(2,3)]
colorRamp <- colorRampPalette(selected_colours)(111)
colorRamp <- colorRampPalette(selected_colours)(46)
colorRamp_f <- colorRampPalette(selected_colours)(6)


#setwd("//gsamba/dataset/RC_Nanopore/itmp/Analysis/results/amrplusplusResis/EnrWaterRes/Rare")
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA")

# method via: https://www.gerkelab.com/blog/2018/09/import-directory-csv-purrr-readr/ #
## Import multiple gene tsv files ##

tsv_files <- fs::dir_ls("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA", regexp = "\\_gene.tsv$")
tsv_files
df <- tsv_files %>% map_dfr(read_tsv, col_names = FALSE, .id = "sample") %>% as_tibble()

## No column names in the original files. Add sample as id and add 3 column names as specified ##
colnames(df) <- c("Sample", "PropReads", "NumUniqGene")
df
## Remove long filenames (we had to import these as sample names) ##
## Edit this based on your file name #
df = df %>%  mutate(Sample = str_remove(Sample, "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA"))
df2 = df %>%  mutate(Sample = str_remove(Sample, ".gene.tsv"))
#df2 = df2 %>%  mutate(Sample = str_remove(Sample, "/rarePE_"))
df2 = df2 %>%  mutate(Sample = str_remove(Sample, ".filt.sorted"))
df2 = df2 %>%  mutate(Sample = str_remove(Sample, "/trim_"))
df2

library(readxl)
## These are the water samples so I don't need to worry about including sample metadata ##
## But if you wanted to include metadata, use a similar code as below and then merge dataframes by sample ##
## Add Site and Visit metadata from lenDat dataframe ##

#lenDat <- read_csv("Filepath/Metadata.csv", col_names = TRUE)
#df3 = df2 %>% left_join(select(lenDat, Sample, Site, Visit), by = c("Sample" = "Sample"))
lenDat <- read_xlsx("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA/Metadata_Rare.xlsx", col_names = TRUE)
lenDat
df3<-merge(df2,lenDat, by = "Sample", all.x = T, all.y = T)
#df3 = df2 %>% left_join(select(lenDat, SampleType, Visit, HerdID), by = c("Sample" = "Sample"))
df3

#If you need to see which files didn't have enough reads to be considered for the graphs.
na_rows <- df3[!complete.cases(df3), ]

## Now let's plot the data using ggplot ##

ggplot(df3, aes(x = PropReads, y = NumUniqGene, color = Sample)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_line(size=1.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  scale_color_manual(values=colorRamp) + 
  facet_wrap(~SampleType)


ggplot(df3, aes(x = PropReads, y = NumUniqGene)) +
  theme_bw() +
  geom_line(aes(color = FarmID), size = 1.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  facet_wrap(~ SampleType) +
  scale_color_manual(values = colorRamp_f)

ggplot(df3, aes(x = PropReads, y = NumUniqGene)) +
  theme_bw() +
  geom_line(aes(color = FarmID), size = 1.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  facet_wrap(~ SampleType) +
  scale_color_brewer(palette = "Set2")

                      

ggplot(df3, aes(x = PropReads, y = NumUniqGene, color = Sample)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_line(size=1.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  scale_color_manual(values=c(colorRamp))


## Now lets look closer at some of the samples with very high genes ##

## define the interesting samples and their colours ##
#cols <- c("H05" = "darkblue", "H10" = "darkred")

df3 <- df3 %>%
  mutate(Category = case_when(
    SampleType == "Bobby F" ~ "Non-replacement calves' faeces",
    SampleType == "Keeper F" ~ "Replacement calves' faeces",
    SampleType == "Replacement M" ~ "Replacement milks",
    SampleType == "Waste M" ~ "Waste milks",
    TRUE ~ SampleType  # Keep all other values unchanged
  ))

# Update the cols object to use new names
cols <- c("Farm3 time 2" = "darkred", "Farm6 time 2" = "darkblue")  # Use new names

# Or if you have an existing cols object, rename it:
names(cols)[names(cols) == "H05"] <- "Farm3 time 2"
names(cols)[names(cols) == "H10"] <- "Farm6 time 2"

# Change specific sample names
df3<- df3 %>%
  mutate(
    Sample = case_when(
      Sample == "H05" ~ "Farm3 time 2",
      Sample == "H10" ~ "Farm6 time 2",
      TRUE ~ Sample  # Keep all other samples unchanged
    )
  )

df3 %>%
  mutate(interesting = case_when(Sample %in% names(cols) ~ Sample,
                                 TRUE ~ NA_character_)) %>%
  ggplot(aes(x = PropReads, y = NumUniqGene, color = interesting, group = Sample)) +
  theme_bw() +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=11)) +
  guides(linetype = guide_legend(override.aes = list(size = 4))) +
  geom_line(size=1.1) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  facet_grid(. ~ Category) +
  scale_color_manual(values= cols,
                     breaks = names(cols),
                     name="Samples",
                     na.value = "gray")
ggsave("../Rarefaction_Oct21.pdf", dpi = "print", units = "cm", width = 26, height = 20)


#setwd("//gsamba/dataset/RC_Nanopore/itmp/Analysis/results/amrplusplusResis/EnrWaterRes/Rare")
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA/")

# method via: https://www.gerkelab.com/blog/2018/09/import-directory-csv-purrr-readr/ #
## Import multiple gene tsv files ##

tsv_files_all <- fs::dir_ls("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA/", regexp = "\\_gene.tsv$")
tsv_files_all
df_all <- tsv_files_all %>% map_dfr(read_tsv, col_names = FALSE, .id = "Sample") %>% as_tibble()

## No column names in the original files. Add sample as id and add 3 column names as specified ##
colnames(df_all) <- c("Sample", "PropReads", "NumUniqGene")
df_all
## Remove long filenames (we had to import these as sample names) ##
## Edit this based on your file name #
df_all = df_all %>%  mutate(Sample = str_remove(Sample, "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA/"))
df2_all = df_all %>%  mutate(Sample = str_remove(Sample, ".gene.tsv"))
df2_all = df2_all %>%  mutate(Sample = str_remove(Sample, "/rarePE_"))
df2_all = df2_all %>%  mutate(Sample = str_remove(Sample, ".filt.sorted"))
df2_all = df2_all %>%  mutate(Sample = str_remove(Sample, "trim_"))
df2_all

library(readxl)
## These are the water samples so I don't need to worry about including sample metadata ##
## But if you wanted to include metadata, use a similar code as below and then merge dataframes by sample ##
## Add Site and Visit metadata from lenDat dataframe ##

#lenDat <- read_csv("Filepath/Metadata.csv", col_names = TRUE)
#df3 = df2 %>% left_join(select(lenDat, Sample, Site, Visit), by = c("Sample" = "Sample"))
lenDat_all <- read_xlsx("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 1 - Cross Sectional Study/R/BWA/Metadata_Rare.xlsx", col_names = TRUE)
lenDat_all
df3_all<-merge(df2_all,lenDat_all, by = "Sample", all.x = T, all.y = T)
#df3 = df2 %>% left_join(select(lenDat, SampleType, Visit, HerdID), by = c("Sample" = "Sample"))
df3_all

## Now let's plot the data using ggplot ##


ggplot(df3_all, aes(x = PropReads, y = NumUniqGene)) +
  theme_bw() +
  geom_line(aes(color = SampleType), size = 1.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  facet_wrap(~ SampleType) +
  scale_color_brewer(palette = "Set2")

ggplot(df3_all, aes(x = PropReads, y = NumUniqGene, group = interaction(FarmID, Sample))) +
  theme_bw() +
  geom_line(aes(color = FarmID), size = 1.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Proportion of reads sampled") + ylab("Number of genes identified") +
  facet_wrap(~ SampleType) +
  scale_color_brewer(palette = "Set2")






