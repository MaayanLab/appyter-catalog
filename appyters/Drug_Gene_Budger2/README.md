# Dr. Gene Budger (DGB) 2

The Dr. Gene Budger 2 (DGB2) Appyter takes a single human gene as input, and returns ranked up- and down-regulating drugs from three Connectivity Mapping resources that were shown to maximally increase or decrease the mRNA expression of the gene in human cell lines. The three Connectivity Mapping resources are:

- [Ginkgo GDPx1 and GDPx2 datasets](https://huggingface.co/ginkgo-datapoints)

- [Novartis DRUG-seq U2OS MoABox dataset](https://zenodo.org/records/14291446)

- [LINCS L1000 Chemical Perturbation dataset](https://maayanlab.cloud/sigcom-lincs/#/Download)

In addition to producing tables of ranked up- and down-regulating drugs of the input gene, the notebook creates volcano plot visualizations and UpSet plots that identify overlap in regulators across datasets. 