install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("VGAM")
install.packages("igraph")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release")

