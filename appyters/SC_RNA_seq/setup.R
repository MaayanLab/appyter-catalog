install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("VGAM")
install.packages("igraph")

install.packages("usethis")
install.packages("gh")

install.packages("https://cran.r-project.org/src/contrib/Archive/XML/XML_3.98-1.20.tar.gz")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release")

install.packages("hdf5r")
devtools::install_github("cran/SDMTools")
install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-74.tar.gz")
install.packages("Hmisc")
install.packages("https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-6.tar.gz")

BiocManager::install("multtest")
install.packages("mutoss")
install.packages("metap")
# install.packages("https://cran.r-project.org/src/contrib/Archive/mutoss/mutoss_0.1-10.tar.gz")
# install.packages("https://cran.r-project.org/src/contrib/Archive/metap/metap_1.2.tar.gz")
source("https://z.umn.edu/archived-seurat")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
devtools::install_github("BaderLab/Tempora")