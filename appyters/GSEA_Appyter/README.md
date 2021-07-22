# GSEA Appyter

The GSEA Appyter performs gene set enrichment analysis and visualizes the results. 

## Inputs
Users must submit gene expression data and a gene set library for the appyter to run. For the gene expression submissions, users have the option of choosing to input either a pre-ranked gene list, or a gene expression dataset with a matching phenotype labels file and a choice of ranking method. As for submitting a gene set library, users can either choose one of the Enrichr libraries, or upload their own. 
Follow this <a href="https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats">link</a> to learn more about the proper file formats.

After determining the data submissions, users can customize their experience further by choosing what enrichment statistic criteria they want their gene sets sorted by, and how many of these top results will be displayed by the appyter.

## Results
A table will display the enrichment statistics for the top gene sets chosen by the input criteria (with the full results available for download as a CSV file), and an interactive Enrichment Plot and Hit Indices Plot pair will also be generated for these sets (available for download as a PNG file).
