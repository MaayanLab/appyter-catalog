# hTFtarget ETL Appyter
### Authors
Ido Diamant - Bioinformatics Software Engineer

Ma’ayan Lab, Mount Sinai Center for Bioinformatics, Department of Pharmacological Sciences  
Icahn School of Medicine at Mount Sinai, New York, NY 10029 USA
### hTFtarget 2022 Dataset
**Genes: 24455**  
**Terms: 1710**  
**Data Source:** http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/

[hTFtarget](http://bioinfo.life.hust.edu.cn/hTFtarget) is a database of human transcription factors. It provides tools for the visualization, interpretation, and analysis of pathway knowledge.

This appyter takes data from the hTFtarget human transcription factor database and outputs files that are usable for Machine Learning and other applications. It processes the [TF-Target-information.txt](http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt) file downloaded on 09-22-2022.
  
The Appyter uses the NCBI database to map the gene names to a set of approved gene symbols so that synonymous genes are mapped to the same symbol.

The Appyter creates gene and attribute similarity matrices, which contain the Jaccard Index between any two genes or attributes.
    
The following output files are made available for download:  
* A binary matrix
* Gene list
* Attribute list
* A gene set library: for each attribute (pathway), a list of genes that are associated with the attribute
* An attribute set library: for each gene, a list of attributes (TFs and tissues) that are associated with each gene
* Gene-gene similarity matrix
* Attribute-attribute similarity matrix
* Gene-attribute edge list: a list of gene-attribute pairs and the strength of each 
association
* Serialized data for Knowledge Graph ingestion: a list of gene and TF:tissue nodes, and gene &rarr; TF:Tissue edges  
  
A ZIP archive containing these files is provided at the bottom of the report.