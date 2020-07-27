# L1000FWD-KD

The L1000KD is a appyter that allows the user to explore 24,000 gene expression signatures resulting from shRNA gene knockdown. 

When the user enters the appyter, they have the choice to input the name of a disease, small molecule or gene, or they are allowed to input their own gene signature; if the user chooses the former, a gene signature of the input is retrieved. After this selection, the inputted signature is compared to all 24,000 signatures, and the user is shown the most similar and opposite gene knockdown signatures, and are given a similarity score, z-score, and p-value for each relevant signature. Next, enrichment analysis is performed on a gene set of the knockdown genes in the 50 most similar signatures. 

Users are also presented with an interactive UMAP dimensional reduction visualization of all of the genes. Users can choose to color this by the name of the knockdown gene, chromosome the gene is on, KEGG pathway, GO pathway, cell line, etc. Users are encouraged to zoom in on specific clusters and examine the contents for potential gene-gene relationships. 

Overall, using a dataset of gene knockdown signatures allows users to investigate the unique functions of individual genes. The L1000KD appyter presents a novel platform to gain insight into gene-to-gene, small molecule-gene, and disease-gene relationships and has promise for generating real-world hypotheses in minutes.