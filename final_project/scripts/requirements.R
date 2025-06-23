# CRAN packages
install.packages(c(
  "tidyverse=2.0.0",
  "ggplot2=3.4.2",
  "here=1.0.1",
  "readxl=1.4.2",
  "scico=1.3.1",
  "data.table=1.14.8",
  "loomR=0.2.1.9000",
  "ggpubr=0.6.0",
  "pheatmap=1.0.12",
  "VennDiagram=1.7.3",
  "clusterProfiler=4.8.1",
  "enrichplot=1.20.0",
  "readr=2.1.4"
), dependencies = TRUE)

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "org.Hs.eg.db=3.17.0"
), version = "3.17")