# Drugmonizome ETL: PharmGKB

[PharmGKB](https://www.pharmgkb.org/) is a pharmacogenomics knowledge resource that encompasses clinical information including clinical guidelines and drug labels, potentially clinically actionable gene-drug associations, and genotype-phenotype relationships.

This appyter takes data from the PharmGKB and outputs files that are usable for Drugmonizome. Drug-gene and drug-variant association data is processed in order to construct a binary matrix with small molecules as rows and genes as column attributes. From this matrix, drug and attribute set libraries are constructed.

Additionally, gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes, are created.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecules
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix