# Harmonizome ETL: The Human Protein Atlas (RNA-Seq)

[The Human Protein Atlas](https://www.proteinatlas.org/) (THPA) is a program with the aim to map all the human proteins in cells, tissues and organs using an integration of various omics technologies. It consists of six parts, each focusing on a different aspect of the proteonome: the Tissue, Cell, Pathology, Blood, Brain, and Metbolic Atlases.

This appyter takes data from RNA-Seq data of human tissues or blood cells and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct an expression matrix with Ensembl gene IDs as rows and tissues or cell lines as column attributes. It then draws from the current NCBI database to map the Ensembl gene IDs to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

From here, it merges duplicate genes and attributes by averaging the expression value and removes genes and attributes that have more than 95% missing data, imputing any remaining missing data for each gene by taking its average expression. It performs normalization by log2 transformation and then quantile normalization for the attributes, and z-scores for the genes.

From this normalized matrix, it creates a matrix of standardized values between -1 and 1. From this standard matrix, it creates a ternary matrix, which stores whether a gene and attribute are negatively, positively, or not correlated. It also creates gene and attribute similarity matrices, which store the cosine distance between any two genes or attributes.

The downloadable file will have the following outputs:
* Unfiltered matrix: the expression matrix before normalization
* Filtered matrix: the normalized matrix
* Gene list
* Attribute list
* Standard matrix
* Ternary matrix
* Up/down gene set library: for each attribute, a list of genes that are positively/negatively correlated, respectively
* Up/down attribute set library: for each gene, a list of attributes that are positively/negatively correlated, respectively
* Gene similarity matrix
* Attribute similarity matrix
* Gene-attribute edge list: a list of gene-attribute pairs and the standard value for each pair 
