# Drugmonizome ETL: RDKit

[RDKit](https://www.rdkit.org/) is an open source chemoinformatics tool used to generate molecular fingerprints.

This appyter converts SMILES string representations of small molecules into chemical fingerprints using RDKit and outputs files that are usable for Drugmonizome. These data are processed in order to construct a binary matrix with small molecules as rows and chemical fingerprints as column attributes. From this matrix, drug and attribute set libraries are constructed.

Additionally, gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes, are created.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecules
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix