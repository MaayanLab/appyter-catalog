# L1000FWD Consensus Appyter

The LINCS Data Portal 3.0 (LDP 3.0) hosts ranked L1000 [1] perturbation signatures from a variety of perturbation types including: drugs and other small molecules, CRISPR knockouts, shRNA knockdowns, and single overexpression. LDP 3.0's RESTful APIs enable querying the signatures programmatically to identify mimickers or reversers for input up and down gene sets. This appyter extends this functionality by enabling analysis for a collection of input signatures to identify consistently reoccuring mimickers and reversers. The appyter takes as input a set of two-sided or one-sided gene sets and constructs a score matrix of mimicking and reversing signatures. From this matrix the appyter computes the consensus. The pipeline also includes (1) Clustergrammer [2] interactive heatmap, and (2) enrichment analysis of the top gene perturbations [3-6] to elucidate the pathways that are being targeted by the consensus perturbagens.

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
1] Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., Natoli, T. E., Lu, X., ... & Golub, T. R. (2017). A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. Cell, 171(6), 1437-1452.

[2] Fernandez, N. F., Gundersen, G. W., Rahman, A., Grimes, M. L., Rikova, K., Hornbeck, P., & Ma’ayan, A. (2017). Clustergrammer, a web-based heatmap visualization and analysis tool for high-dimensional biological data. Scientific data, 4(1), 1-12.

[3] Chen, E. Y., Tan, C. M., Kou, Y., Duan, Q., Wang, Z., Meirelles, G. V., ... & Ma’ayan, A. (2013). Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC bioinformatics, 14(1), 1-14.

[4] Kuleshov, Maxim V., et al. "Enrichr: a comprehensive gene set enrichment analysis web server 2016 update." Nucleic acids research 44.W1 (2016): W90-W97.

[5] Xie, Z., Bailey, A., Kuleshov, M. V., Clarke, D. J., Evangelista, J. E., Jenkins, S. L., ... & Ma'ayan, A. (2021). Gene set knowledge discovery with Enrichr. Current protocols, 1(3), e90.

[6] Kropiwnicki, E., Evangelista, J. E., Stein, D. J., Clarke, D. J., Lachmann, A., Kuleshov, M. V., ... & Ma’ayan, A. (2021). Drugmonizome and Drugmonizome-ML: integration and abstraction of small molecule attributes for drug enrichment analysis and machine learning. Database, 2021.

[7] Maglott, D., Ostell, J., Pruitt, K. D., & Tatusova, T. (2005). Entrez Gene: gene-centered information at NCBI. Nucleic acids research, 33(suppl_1), D54-D58.

[8] Brown, G. R., Hem, V., Katz, K. S., Ovetsky, M., Wallin, C., Ermolaeva, O., ... & Murphy, T. D. (2015). Gene: a gene-centered information resource at NCBI. Nucleic acids research, 43(D1), D36-D42.