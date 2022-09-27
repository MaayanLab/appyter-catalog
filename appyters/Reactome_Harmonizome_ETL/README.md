# Harmonizome ETL: Reactome

[Reactome](https://reactome.org/) is a database of manually curated pathways. It provides tools for the visualization, interpretation, and analysis of pathway knowledge.

This appyter takes data from the Reactome Pathways Gene Set and outputs files that are usable for Machine Learning and other applications. It processes the [ReactomePathways.gmt.zip](https://reactome.org/download/current/ReactomePathways.gmt.zip) file.
  
The Appyter uses the NCBI database to map the gene names to a set of approved gene symbols so that synonymous genes are mapped to the same symbol.

The Appyter creates gene and attribute similarity matrices, which contain the Jaccard Index between any two gene sets or attribute sets.
    
The following output files are made available for download:
* A binary matrix
* Gene list
* Attribute list
* A gene set library: for each attribute (pathway), a list of genes that are associated with the attribute
* An attribute set library: for each gene, a list of attributes (pathways) that are associated with each gene
* Gene-gene similarity matrix
* Attribute-attribute similarity matrix
* Gene-attribute edge list: a list of gene-attribute pairs and the strength of each 
association
* Serialized data for Knowledge Graph ingestion: a list of gene and pathway nodes, and gene &rarr; pathway edges

A ZIP archive containing these files is provided at the bottom of the report.
