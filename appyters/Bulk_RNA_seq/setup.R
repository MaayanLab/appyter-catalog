if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")

install.packages("statmod")