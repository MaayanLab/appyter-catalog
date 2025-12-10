# DrugRanger

The DrugRanger Appyter takes a single human gene or gene set as input, and returns ranked up- and down-regulating drugs from three Connectivity Mapping resources that were shown to maximally increase or decrease the mRNA expression of the gene(s) in human cell lines. The three Connectivity Mapping resources are:

- [Ginkgo GDPx1 and GDPx2 datasets](https://huggingface.co/ginkgo-datapoints)

- [Novartis DRUG-seq U2OS MoABox dataset](https://zenodo.org/records/14291446)

- [Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M)

In addition to producing tables of ranked up- and down-regulating drugs of the input gene, the notebook creates various visualizations for the single gene and multi-gene analysis, to help users determine the most effective regulators of their input gene(s). 