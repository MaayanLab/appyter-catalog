# Tumor Gene Target Screener

This Appyter is inspired by the work of Bosse, Kristopher R et al [1] which compared neurobastomas vs normal tissue in GTEx [2] to identify a promising candidate immunotherapeutic target.

The goal is to allow rapid screening of targets with the help of normal tissue data from GTEx and GEO data through ARCHS4 [3], as well as single-cell data from Tabula Sapiens and the Human Cell Atlas. The Appyter takes tumor expression data and attempts to rank significantly differentially expressed genes when compared with with either bulk RNA-seq data from GTEx or ARCHS4, or single-cell RNA-seq data from Tabula Sapiens or Human Cell Atlas, across all tissues.

The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of the Office of the Director of the National Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. [2] GTEx Version 8 gene counts was processed to produce gene summary statistics.

ARCHS4 provides access to gene counts from HiSeq 2000, HiSeq 2500 and NextSeq 500 platforms for human and mouse experiments from GEO and SRA. [3] We processed ARCHS4 Version 11 to produce gene summary statistics.

The [Tabula Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org) dataset was created by the The Tabula Sapiens Consortium. [4] We processed the Tabula Sapiens dataset to produce gene summary statistics.

The [Human Cell Atlas](https://data.humancellatlas.org) provides access to single-cell data contributed by the scientific community.  We combined and processed 15 datasets from the Human Cell Atlas to produce gene summary statistics.

Immunotherapeutic candidates must have limited expression in normal tissues to be considered safe targets, so proteomic visualizations of the highly expressed genes in normal tissues may be useful in assessing gene candidacy. Proteomics data were obtained from the [Human Protein Atlas](https://www.proteinatlas.org/about/download) [5] with IHC-based expression profiling, the [Human Proteome Map](https://www.humanproteomemap.org/download.php) [6] with MS-based expression quantification, and a [GTEx proteome project](https://doi.org/10.1016/j.cell.2020.08.036) [7] using TMT MS. 


[1] Bosse, Kristopher R et al. "Identification of GPC2 as an Oncoprotein and Candidate Immunotherapeutic Target in High-Risk Neuroblastoma." Cancer cell vol. 32,3 (2017): 295-309.e12. <https://doi.org/10.1016/j.ccell.2017.08.003>

[2] Lonsdale, John, et al. "The genotype-tissue expression (GTEx) project." Nature genetics 45.6 (2013): 580-585. <https://doi.org/10.1038/ng.265>

[3] Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, Ma'ayan A. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications 9. Article number: 1366 (2018), <https://doi.org/10.1038/s41467-018-03751-6>

[4] The Tabula Sapiens Consortium, Science 376, eabl4896 (2022)

[5] Uhlén M et al. "Tissue-based map of the human proteome." Science (New York, N.Y.) vol. 347,6220 (2015): 1260419. <https://doi.org/10.1126/science.1260419>

[6] Kim, Min-Sik et al. “A draft map of the human proteome.” Nature vol. 509,7502 (2014): 575-81. <https://doi.org/10.1038/nature13302>

[7] Jiang, Lihua et al. “A Quantitative Proteome Map of the Human Body.” Cell vol. 183,1 (2020): 269-283.e19. <https://doi.org/10.1016/j.cell.2020.08.036>