# Harmonizome ETL: COSMIC (Copy Number Variants)

[The Catalogue of Somatic Mutations in Cancer](https://cancer.sanger.ac.uk/cosmic) (COSMIC) is an expert-curated database of somatic mutation information related to human cancers. It combines high-precision data from peer reviewed papers with genome-wide screen data to provide extensive converage of the cancer genomic landscape.

This appyter takes data from the COSMIC Complete Mutation Data and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct a binary matrix with gene names as rows and primary sites as column attributes. It then draws from the current NCBI database to map the gene names to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

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
