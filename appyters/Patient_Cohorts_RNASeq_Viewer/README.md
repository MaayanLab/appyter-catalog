# Patient Cohorts RNA-Seq Viewer

![Preview image](Patient_Cohorts_RNASeq_Viewer/static/main-image.png)

[The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) dataset contains multiomics profiling and clinical data from over 10,000 tumors collected from patients spanning several cancer types. Specifically, TCGA has bulk RNA-sequencing (RNA-Seq) profiling of tumors, which can provide insights into mechanisms and classify tumors by subtype.

By default, this appyter provides analysis and visualization of TCGA datasets. Users can optionally upload their own datasets.

The appyter provides analysis for RNA-Seq TCGA data for cancers with over 150 cases. The report automatically identifies clusters of patient and determines which clinical features and genes are most associated with each cluster.

For the TCGA data, each column in the RNA-Seq dataset corresponds to a row in the clinical dataset; both are referenced by the same identifier (here the case_id as provided by TCGA).

The RNA-Seq data loaded from TCGA is in the form of raw counts mapped to genes with the [htseq-count](https://htseq.readthedocs.io/en/release_0.9.0/count.html) analysis package; the same format should be followed for user-uploaded files. The analysis filters out lowly expressed genes, identifies the most variable genes, normalize the counts, and reduces the dimensionality of the dataset further with PCA and UMAP.

To determine the ideal number of clusters, the analysis tests a range of possible K clusters, and selects the optimal number based on a modified silhouette score that prioritizes more clusters to avoid missing out small clusters.

The appyter also identifies the top genes for each cluster, using these for enrichment analysis and suggestion for drugs and small molecules based on the drugs that mimic or reverse the signatures obtained for each cluster. Such drug suggestions are based on the L1000 dataset, using the L1000FWD API. It should be noted that these are speculative predictions and should not be applied to patients before carefully tested in cell based assays and animal models.
