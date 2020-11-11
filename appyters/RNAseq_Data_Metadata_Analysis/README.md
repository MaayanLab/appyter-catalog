# RNA-seq and Metadata Analysis Appyter

This notebook template provides a pipeline for the visualization and analysis of RNA-seq gene read counts and associated metadata. The analysis process requires a rich metadata file with multiple categorical columns, but is largely generalizable. If using GEO data with only one metadata column containing tags, please use the [Bulk RNA-seq Analysis Appyter](https://appyters.maayanlab.cloud/Bulk_RNA_seq/). If interested in GTEx data analyses, please use the GTEx Tissue-Specific RNA-seq Analysis Appyter. 

The RNA-seq data is normalized before undergoing Principle Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP). The data is then clustered using k-means clustering, with `k` determined by calculating a silhouette score. Clusters are visualized using the [React-Scatter-Board](https://github.com/MaayanLab/react-scatter-board) package. 

For each cluster, the most up-regulated and down-regulated genes are also identified and used for enrichment analysis via the [Enrichr](https://maayanlab.cloud/Enrichr/) API. Enrichment results are visualized with the [React-GSEA](https://github.com/MaayanLab/react-GSEA/tree/simplified) package.

Similar and opposite drug/small molecule perturbation signatures from the L1000 database can be queried for each cluster using the [L1000FWD](https://maayanlab.cloud/L1000FWD/) API. This step is optional, but may be useful for disease-specific RNA-seq data, for which the goal is to identify potential treatments. 

Note: If using GTEx data or other healthy tissue sample data for which querying drug signatures is not relevant, please use the GTEx Tissue-Specific RNA-seq Analysis Appyter instead. If using GEO data, please use the Bulk RNA-seq Analysis Appyter.