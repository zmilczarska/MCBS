# fig_3b.R
# top 25 clusters for scz2022

library(data.table)
library(tidyverse)
library(ggplot2)
library(scico)
library(readxl) 
# read in SCZ cluster level results
df <- read_xlsx("../data/Supplementary_Datasets.xlsx", sheet = "Supplementary_Data_8") |> 
  filter(Phenotype == "scz2022") |> 
  slice_min(P,  n = 25)


p <- ggplot(
  df,
  aes(
    x = reorder(Cluster, -P.fdr),
    y = -log10(P),
    fill = -log10(P)
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("") + ylab(expression('-log'[10]*'(P)')) +
  labs(title="Top Clusters (scz2022)") +
  theme_light() + theme(axis.title.x=element_text(hjust=0)) +
  scico::scale_fill_scico(palette = "lajolla",direction = -1, midpoint = 5) +
  guides(fill="none")

ggsave(p,
  file = "../figures/3b_topClusters_scz.pdf",
  width = 5, height = 6
)

# read in neurotism cluster level results
df_1 <- read_xlsx("../data/Supplementary_Datasets.xlsx", sheet = "Supplementary_Data_8") |> 
  filter(Phenotype == "neuroticism") |> 
  slice_min(P,  n = 25)

p_1 <- ggplot(
  df_1,
  aes(
    x = reorder(Cluster, -P.fdr),
    y = -log10(P),
    fill = -log10(P)
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("") + ylab(expression('-log'[10]*'(P)')) +
  labs(title="Top Clusters (neuroticism)") +
  theme_light() + theme(axis.title.x=element_text(hjust=0)) +
  scico::scale_fill_scico(palette = "lajolla",direction = -1, midpoint = 5) +
  guides(fill="none")

ggsave(p_1,
       file = "../figures/3b_topClusters_neurotism.pdf",
       width = 5, height = 6
)
# read in educattional attaint. cluster level results
df_2 <- read_xlsx("../data/Supplementary_Datasets.xlsx", sheet = "Supplementary_Data_8") |> 
  filter(Phenotype == "educational_attainment") |> 
  slice_min(P,  n = 25)

p_2 <- ggplot(
  df_2,
  aes(
    x = reorder(Cluster, -P.fdr),
    y = -log10(P),
    fill = -log10(P)
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("") + ylab(expression('-log'[10]*'(P)')) +
  labs(title="Top Clusters (edu)") +
  theme_light() + theme(axis.title.x=element_text(hjust=0)) +
  scico::scale_fill_scico(palette = "lajolla",direction = -1, midpoint = 5) +
  guides(fill="none")

ggsave(p_2,
       file = "../figures/3b_topClusters_edu.pdf",
       width = 5, height = 6
)
