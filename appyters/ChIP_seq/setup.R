if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")

BiocManager::install("ChIPseeker")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")