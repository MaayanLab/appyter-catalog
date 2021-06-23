# ChEA3 Consensus Kinases

This Appyter takes as input a collection of gene sets and returns the top consensus transcription factorss associated with them. The Appyter performs transcription factors enrichment analysis for all uploaded gene sets using [ChEA3](https://maayanlab.cloud/ChEA3/). The top consensus transcription factorss are determined using mean rank values. The enriched transcription factorss are then visualized using downloadable stacked bar plots and static heatmaps.

## File Format
The appyter takes up and down gene sets gmt files with the following format:
```
signature 1\t\tGENE 1\tGENE 2\tGENE 3
signature 2\t\tGENE 3\tGENE 4\tGENE 5
signature 3\t\tGENE 2\tGENE 5\tGENE 6
signature 4\t\tGENE 1\tGENE 2\tGENE 6
```

Signature names should be _the same_ on both up and down gmt files. Below are the sample up and down gene sets:

* [Sample Gene Set](https://appyters.maayanlab.cloud/storage/ChEA3Consensus/sars_cov_2_down_genesets.gmt)
* [Sample Meta File](https://appyters.maayanlab.cloud/storage/ChEA3Consensus/sars_cov_2_meta.tsv)
