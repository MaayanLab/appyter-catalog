if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz')
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")

install.packages("statmod")