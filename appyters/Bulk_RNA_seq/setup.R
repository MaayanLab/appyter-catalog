if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = '3.9', ask = FALSE)

BiocManager::install("limma")

install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz')
BiocManager::install("edgeR")

BiocManager::install("DESeq2")

install.packages("statmod")

# verify
require(limma)
require(edgeR)
require(DESeq2)
