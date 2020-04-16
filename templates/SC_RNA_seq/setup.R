if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("monocle")
devtools::install_github('LTLA/SingleR')
BiocManager::install("scater")
BiocManager::install("scRNAseq")


install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("VGAM")
install.packages("igraph")
