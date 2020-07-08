# Harmonizome ETL: Public Health Genomics and Precision Health Knowledge Base

[the Public Health Genomics and Precision Health Knowledge Base](https://phgkb.cdc.gov/PHGKB/) (PHGKB) a database of published scientific literature, CDC resources, and other materials about genomics and precision health discoveries.

This appyter takes data from the Genopedia and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct a binary matrix with gene names as rows and associated diseases as column attributes. It then draws from the current NCBI database to map the gene names to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

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
