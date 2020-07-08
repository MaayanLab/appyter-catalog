# Harmonizome ETL: The Human Protein Atlas (Immunihistochemistry)

[The Human Protein Atlas](https://www.proteinatlas.org/) (THPA) is a program with the aim to map all the human proteins in cells, tissues and organs using an integration of various omics technologies. It consists of six parts, each focusing on a different aspect of the proteonome: the Tissue, Cell, Pathology, Blood, Brain, and Metbolic Atlases.

This appyter takes data from expression profiles based on immunohistochemistry and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct a binary matrix with Ensembl gene IDs as rows and tissue and cell types as column attributes. It then draws from the current NCBI database to map the Ensembl gene IDs to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

From here, it creates gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes. 

The downloadable file will have the following outputs:
* Binary matrix: the expression matrix with gene symbols
* Filtered matrix: the normalized matrix
* Gene list
* Attribute list 
* Up gene set library: for each attribute, a list of genes that are correlated
* Up attribute set library: for each gene, a list of attributes that are correlated
* Gene similarity matrix
* Attribute similarity matrix
* Gene-attribute edge list: a list of gene-attribute pairs and the expression for each pair 
