This appyter is the first part of the two-part ENKEFALOS analysis pipeline. It takes a list of differentially expressed genes (recommended to be from a neural tissue/cell sample) and helps identify genes that are significantly associated with neural electrophysiological and/or morphological measures using data derived from a study done by the Pavlidis Lab at the University of British Columbia. You can find the description of the data and how it was derived [here](https://github.com/PavlidisLab/transcriptomic_correlates). The second part of the ENKEFALOS appyter can be utilized for more downstream, single-gene analysis based on your results from this appyter's results.

# Framework

<p align="center">
  <img width="615" alt="image" src="https://github.com/KrishU27/Enkefalos/assets/132734331/c378127d-4168-43cd-8907-46f2e0a65e3f">
</p>

There are several sections in this appyter, for which we have a brief overview below. If you would like a more comprehensive guide to how to use ENKEFALOS, please refer to our user guide [here](https://docs.google.com/document/d/15h8A65FygTK2_KLA_-6R8u8clJSRnF7yHLFAvN7BN6Y/edit?tab=t.0#heading=h.ola4n01ccsle).

Appyter 1:
- Takes in your genes of interest (GOI) as well as a FDR threshold for analyses.
- Displays genes from your GOI that had significant correlations with electrophysiological/morphological measures.
- Prints out the number of enriched genes as well as what the genes with significant correlations are.
- Calls StringDB to create a gene interactome of your enriched genes. Will tabulate the number of interactions each gene has.

# Reference
1. [Bomkamp C, Tripathy SJ, Bengtsson Gonzales C, Hjerling-Leffler J, Craig AM, et al. (2019) Transcriptomic correlates of electrophysiological and morphological diversity within and across excitatory and inhibitory neuron classes. PLOS Computational Biology 15(6): e1007113.](https://journals.plos.org/ploscompbiol/article/citation?id=10.1371/journal.pcbi.1007113)
