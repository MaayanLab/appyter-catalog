# Drug Gene Budger (DGB) 2

This appyter takes a single gene as input and identifies up-regulating and down-regulating drugs from three connectivity mapping resources.

- [Ginkgo GDPx1 and GDPx2 datasets](https://huggingface.co/ginkgo-datapoints)

- [Novartis DRUG-seq U2OS MoABox dataset](https://zenodo.org/records/14291446)

- [LINCS L1000 Chemical Perturbation dataset](https://maayanlab.cloud/sigcom-lincs/#/Download)

In addition to producing tables of ranked up- and down-regulating drugs of the input gene, the notebook creates volcano plot visualizations and UpSet plots that identify overlap in regulators across datasets. 