# STEAP

**S**ingle cell **T**ype **E**nrichment **A**nalysis for **P**henotypes (**STEAP**) uses scRNA-seq data and GWAS summary statistics to determine which cell-types are enriched in the GWAS phenotype. It is an extension to [CELLECT](https://github.com/perslab/CELLECT) and uses [S-LDSC](https://github.com/bulik/ldsc) ([Finucane et al., 2015](https://www.nature.com/articles/ng.3404)), [MAGMA](https://ctg.cncr.nl/software/magma) [(de Leeuw et al., 2015)](https://doi.org/10.1371/journal.pcbi.1004219) and [H-MAGMA](https://github.com/thewonlab/H-MAGMA) [(Sey et al., 2020)](https://doi.org/10.1038/s41593-020-0603-0) for enrichment analysis.

STEAP runs multiple post-processing steps on top of CELLECT output files:

- Gene Set Enrichment Analysis (GSEA)
- Cell-Type Correlation
- Expression Specificity (ES) Gene Correlation

![pipeline](https://raw.githubusercontent.com/erwinerdem/STEAP/master/pipeline.png)


This appyter notebook requires the `prioritization.csv` files from the cell type enrichment analysis. To get these for your phenotype of interest you must first run the [STEAP pipeline](https://github.com/erwinerdem/STEAP).
