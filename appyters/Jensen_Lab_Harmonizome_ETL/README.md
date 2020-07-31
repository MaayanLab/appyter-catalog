# Harmonizome ETL: Jensen Lab

[The Jensen Lab Databases](https://compartments.jensenlab.org/Search) are a set of web resources that integrates evidence on tissue expression, disease-gene interaction, and protein subcellular localization from manually curated literature, high-througput screens, automatic text mining, and other methods. They map the evidence to common protein identifiers and Gene Ontology (GO) terms.

This appyter takes data from any of the databases and outputs files that are usable for the Harmonizome. It pre-processes the raw data  in order to construct a binary matrix with gene names as rows and compartments, diseases, or tissues as column attributes. It then draws from the current NCBI database to map the gene names to a set of approved gene symbols, so that synonymous genes are mapped to the same symbol. 

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
