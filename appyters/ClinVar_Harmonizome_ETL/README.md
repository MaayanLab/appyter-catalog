# Harmonizome ETL: ClinVar

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a public archive of reports of relationships between human variations and phenotypes. It provides information about allelic variations regarding their clinical significance and other supporting data.

This appyter takes data from a summary of allelic variants and outputs files that are usable for the Harmonizome. It pre-processes the raw data by selecting for allelic variants with clinical significance, and which have been reviewed in order to construct a binary matrix with gene names as rows and phenotypes as column attributes. It then draws from the current NCBI database to map the gene names to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

From here, it creates gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes. 

The downloadable file will have the following outputs:
* Binary matrix: the expression matrix with gene symbols
* Gene list
* Attribute list 
* Up gene set library: for each attribute, a list of genes that are correlated
* Up attribute set library: for each gene, a list of attributes that are correlated
* Gene similarity matrix
* Attribute similarity matrix
* Gene-attribute edge list: a list of gene-attribute pairs and the expression for each pair 
