# tcgaEnrichrViewer

[The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) is a massive database containing laboratory analysis and clinical data for over 10,000 cancer patients spanning several hundred cancer subtypes. Specifically, TCGA offers extensive bulk RNA-sequencing (RNA-Seq) data, which can provide insight into the particular genes implicated in specific cancer subtypes.

By default, this appyter uses TCGA data. Users can optionally also upload their own datasets.

The appyter analyzes RNA-Seq data for cancers with over 150 cases in TCGA (or user-provided cohorts), producing clusters of patient subtypes and determining which clinical features and genes are most associated with each cluster.

For TCGA data, each column in the RNA-Seq dataset corresponds to a row in the clinical dataset; both are referenced by the same identifier (here the `case_id` as provided by TCGA).

The RNA-Seq data loaded from TCGA is in the form of raw counts mapped to genes with the [htseq-count](https://htseq.readthedocs.io/en/release_0.9.0/count.html) analysis package; the same format should be obeyed for user-uploaded files. We then filter for the most variable genes, normalize those counts, and reduce the dimensionality of the dataset further with PCA and UMAP.

To determine the ideal number of clusters, we probe a range and select the number based on a modified silhouette score that prioritizes more clusters (so we do not miss out on small clusters).

We identify the top genes for each cluster, using these for Enrichment analysis and treatment suggestion based on the drugs used for perturbation to produce signatures in the L1000 dataset that are most opposite to each cluster.
