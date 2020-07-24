# Drugmonizome ETL: ATC Codes

[DrugBank](https://www.drugbank.ca/) is a bioinformatics and chemoinformatics resource that catalogs extensive information regarding drug and experimental small molecule attributes and metadata.

This appyter extracts data from the DrugBank database and matches fourth level ATC Codes to all unique chemical entities.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecule
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix