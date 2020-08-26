# Patient Cohorts RNA-Seq Viewer

The Patient Cohorts RNA-Seq Viewer Appyter is a web application that provides a customizable interface for processing, visualizing, and analyzing RNA-sequencing (RNA-Seq) data from patient cohorts. The goal of the Appyter is to provide comprehensive analysis of patient cohorts by considering the RNA-seq profiling of the patient samples together with information about their clinical parameters.

The Appyter automatically identifies clusters of patients based on their RNA-seq profiles and associates clinical metadata with each cluster. The Patient Cohorts RNA-Seq Viewer Appyter is preloaded with data collected by The Cancer Genome Atlas (TCGA) but can accommodate user-uploaded datasets.

A standard RNA-seq preprocessing pipeline, including normalization and dimensionality reduction, is implemented. Clusters of patients and their associated differential gene expression profiles are fed into a series of downstream analyses. These include survival analysis, enrichment analysis, and small molecule and drug prioritization based on the L1000 dataset. These analyses provide insights into each cluster’s unique genomic, transcriptomic, and clinical features.

The output of the Appyter is a downloadable Jupyter notebook with code blocks, descriptive markdown, interactive and static figures, and tables. Figure and table legends are provided so users can export the results directly into their publications. Of note, various parameters and methods used for preprocessing and subsequent analyses can be adjusted by the user in the input form before the notebook is executed. The Appyter enables researchers with no programming background to perform complex analyses to uncover patterns embedded in their RNA-seq patient cohort datasets.
