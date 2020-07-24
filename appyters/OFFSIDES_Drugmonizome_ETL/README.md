# Drugmonizome ETL: OFFSIDES

[OFFSIDES](https://purl.stanford.edu/zq918jm7358) is a resource of predicted drug-side effect interactions that were extracted from adverse event reporting systems using an adaptive data-driven approach for correcting cases for which covariates were unknown or unmeasured and combined this approach with existing methods to improve analyses of drug effects.

This appyter extracts data from OFFSIDES and matches small molecules to predicted side-effects based on a p-value threshold specified by the user, and outputs files that are usable for Drugmonizome. Side effect interaction data is processed in order to construct a binary matrix with small molecules as rows and genes as column attributes. From this matrix, drug and attribute set libraries are constructed.

Additionally, gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes, are created.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecule
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix