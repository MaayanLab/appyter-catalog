# Drugmonizome ETL: CREEDS

[CREEDS](https://amp.pharm.mssm.edu/CREEDS/) is a crowdsourcing resource for the curation and reanalysis of gene expression profiles from GEO and includes drug perturbation signatures.

This appyter takes drug-induced gene expression data from CREEDS and outputs files that are usable for Drugmonizome. These data are processed in order to construct upregulated and downregulated binary matrices with small molecules as rows and genes as column attributes. From these matrices, up and down drug and attribute set libraries are constructed.

Additionally, gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes, are created.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecules
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix