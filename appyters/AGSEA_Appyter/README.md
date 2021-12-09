# AGSEA Appyter

The Augmented Gene Set Enrichment Analysis (AGSEA) Appyter performs and visualizes standard gene set enrichment analysis (GSEA), as well as augmented GSEA.

## Inputs
Users must submit gene expression data and a gene set library for the appyter to run. For the gene expression submissions, users have the option of choosing to input either a pre-ranked gene list, or a gene expression dataset with a matching phenotype labels file and a choice of ranking method. As for submitting a gene set library, users can either choose one of the Enrichr libraries, or upload their own. 
Follow this [link](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) to learn more about the proper file formats.

After determining the data submissions, users can customize their experience further by choosing what enrichment statistic criteria they want their gene sets sorted by, how many of these top results will be displayed by the appyter, and whether they want to apply augmented GSEA.

## Results
AGSEA outputs a table containing the enrichment results statistics for the top gene sets from the chosen gene set library selected by the user. The full results are available for download as a CSV file, and an interactive GSEA plot is also generated. The AGSEA Appyter also highlights genes that are co-expressed with each annotated gene set and are highly ranked within the input signature. These genes are typically understudied genes that should be further explored as relevant to the biological process under investigation. 
