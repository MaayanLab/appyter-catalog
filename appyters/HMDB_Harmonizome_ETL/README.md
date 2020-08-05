# Harmonizome ETL: Human Metabolome Database

[The Human Metabolome Database](https://hmdb.ca/) (HMDB) is a database of small molecule metabolites in the human body. It contains chemical, clinical, and biochemistry data. It has applications in metabolomics, clinical chemistry, and biomarker discovery.

This appyter takes data from the complete metabolite dataset and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct a binary matrix with gene names as rows and associated proteins as column attributes. It then draws from the current NCBI database to map the gene names to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

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
