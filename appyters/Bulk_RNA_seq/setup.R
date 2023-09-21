install.packages("R.utils")
install.packages("RCurl")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = '3.16', ask = FALSE)

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
install.packages("statmod")

# verify
require(limma)
require(edgeR)
require(DESeq2)
