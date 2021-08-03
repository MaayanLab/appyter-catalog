# L1000FWD Consensus Drugs

This Appyter takes up and down gene sets and returns the consensus mimicker and reverser perturbagen associated to them. Users can input up and down gene sets to [LINCS Fata Portal 3.0](https://ldp3.cloud/) to get the top mimicker and reverser signatures. It then counts the number of times a perturbagen appear of a query gene set. The appyter only keeps perturbagens that appear in a certain percentage of query gene sets. The p-value for each drug is then computed using a pre-computed empirical distribution using [CREEDS](https://maayanlab.cloud/creeds/#downloads) gene sets as background signatures. The drugs are then visualized using [clustergrammer](https://maayanlab.cloud/clustergrammer/). A publication ready static heatmap is also provided for the users. This appyter also support uploading single direction gene sets (all up or all down gene sets), it then returns mimickers (signatures that ranks the genes in the gene sets similar to its direction) ang reversers (signatures that ranks the genes opposite to the direction of the input gene sets).

## File Format
The appyter takes up and down gene sets gmt files with the following format:
```
signature 1\t\tGENE 1\tGENE 2\tGENE 3
signature 2\t\tGENE 3\tGENE 4\tGENE 5
signature 3\t\tGENE 2\tGENE 5\tGENE 6
signature 4\t\tGENE 1\tGENE 2\tGENE 6
```

Signature names should be _the same_ on both up and down gmt files. Below are the sample up and down gene sets:

* [Up Gene Set](https://appyterbucket.s3.amazonaws.com/sample_data/up_down_signatures/up_diseases)
* [Down Gene Set](https://appyterbucket.s3.amazonaws.com/sample_data/up_down_signatures/down_diseases)

## References
[1] Wang, Z., Monteiro, C. D., Jagodnik, K. M., Fernandez, N. F., Gundersen, G. W., ... & Ma'ayan, A. (2016) Extraction and Analysis of Signatures from the Gene Expression Omnibus by the Crowd. Nature Communications doi: 10.1038/ncomms12846

[2] Fernandez, N. F. et al. Clustergrammer, a web-based heatmap visualization and analysis tool for high-dimensional biological data. Sci. Data 4:170151 doi: 10.1038/sdata.2017.151 (2017).