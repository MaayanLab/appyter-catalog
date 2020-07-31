# Drugmonizome ETL: Geneshot

[Geneshot](https://amp.pharm.mssm.edu/geneshot/) is a resource where users can submit any biomedical search terms to receive prioritized genes that are most relevant to those search terms. Geneshot finds publications that mention both the search terms and genes and outputs associated gene lists for each search term. Additionally, using the associated gene list, further genes can be predicted based on co-occurrence using AutoRIF, GeneRIF, Tagger, and Enrichr or co-expression using ARCHS4.

This appyter takes an input .txt file of small molecule names and outputs files that are usable for Drugmonizome. Associated gene and predicted gene data is processed in order to construct a binary matrix with small molecules as rows and genes as column attributes. From this matrix, drug and attribute set libraries are constructed.

Additionally, gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes, are created.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecules
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix