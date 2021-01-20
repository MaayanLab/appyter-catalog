# Enrichr Consensus Terms

This Appyter takes as input a collection of gene sets and returns the top consensus enrichment terms associated with them. The Appyter performs gene set enrichment analysis for all uploaded gene sets against gene set libraries that the user selects from all the libraries that are available in [Enrichr](https://maayanlab.cloud/Enrichr/). The top consensus terms are determined using p-values returned from the Enrichr enrichment results using the Fisher exact test. The enriched terms are then visualized using downloadable stacked bar plots and static heatmaps.

## File Format
The appyter takes up and down gene sets gmt files with the following format:
```
signature 1\t\tGENE 1\tGENE 2\tGENE 3
signature 2\t\tGENE 3\tGENE 4\tGENE 5
signature 3\t\tGENE 2\tGENE 5\tGENE 6
signature 4\t\tGENE 1\tGENE 2\tGENE 6
```

Signature names should be _the same_ on both up and down gmt files. Below are the sample up and down gene sets:

* [Sample Gene Set](https://appyters.maayanlab.cloud/storage/EnrichrConsensus/sample_input/input.gmt)

## References
[1] Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013;128(14).

[2] Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377.

[3] Fernandez, N. F. et al. Clustergrammer, a web-based heatmap visualization and analysis tool for high-dimensional biological data. Sci. Data 4:170151 doi: 10.1038/sdata.2017.151 (2017).