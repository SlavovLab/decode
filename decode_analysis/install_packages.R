
## all packages required for the amino acid substitution analysis

## from CRAN
install.packages("viridis")
install.packages("readxl")
install.packages("xgboost")
install.packages("ggplot2")
install.packages("dplyr")

## from BioConductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffLogo")

## from github
## https://josephcrispell.github.io/projects/basicplotter (2024-08-24),
## https://github.com/raim/segmenTools (release RAAS_preprint).
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("JosephCrispell/basicPlotteR")
remotes::install_github("raim/segmenTools@RAAS_preprint")
