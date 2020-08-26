# Harmonizome ETL: Bgee

[Bgee](https://bgee.org/) is a database of expression patterns in multiple animal species, produced from multiple data types. It is based on healthy expression data to provide a comparable reference of normal gene expression.

This appyter takes data from Bgee gene expression calls and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct a binary matrix with gene names as rows and a choice of anatomica entity, developmental stage, or both as column attributes. It then draws from the current NCBI database to map the gene names to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

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
