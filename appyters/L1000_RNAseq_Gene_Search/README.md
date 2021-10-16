# Transformed L1000 to RNA-seq Perturbational Signatures Gene Search

This Appyter provides visualizations of the RNA-seq signatures induced by CRISPR knockouts and chemical perturbagens. Signatures are computed from transformed data profiles from the [LINCS L1000 data](). The transformation was performed using a two-step model: 
1. A cycleGAN model was used to first predict the RNA-seq expression of the 978 L1000 landmark genes
2. A fully connected neural network was used to extrapolate the predicted RNA-seq expression of the 978 landmark genes to a full set of 23,164 genes

Signatures were computed using the characteristic direction method [(Clark et al., 2014)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-79), as implemented [here](https://github.com/MaayanLab/maayanlab-bioinformatics/blob/master/maayanlab_bioinformatics/dge/characteristic_direction.py). 