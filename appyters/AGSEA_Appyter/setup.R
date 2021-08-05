install.packages("R.utils")
install.packages("RCurl")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("limma")
BiocManager::install("statmod")
BiocManager::install("edgeR")