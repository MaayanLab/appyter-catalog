# Kinase - Transcription Factor Module Pairwise Analysis

This Appyter predicts shared modules of kinases, inferred from phosphoproteomic data through [KEA3](https://maayanlab.cloud/kea3/) [1], and transcription factors, inferred from transcriptomic data from [ChEA3](https://maayanlab.cloud/chea3/) [2], that are differentially active across two groups of samples. It first ranks the kinases and transcription factors and then scores pairs of kinases and transcription factors on an individual sample level per each group. These scores are binarized based on an additional threshold which only takes into account the top percentage of pair scores and with which clustermaps are created. Clusters or modules of kinases and transcription factors in these groups are then compared and those that appear across multiple of the heatmaps, representing different directions of activation, are retained and their relationship is displayed as a bipartite network.

[1] Kuleshov MV, Xie Z, London ABK, Yang J, Evangelista JE, Lachmann A, Shu I, Torre D, Ma'ayan A. KEA3: improved kinase enrichment analysis via data integration. Nucleic Acids Res. 2021 Jul 2;49(W1):W304-W316. doi: 10.1093/nar/gkab359.

[2] Keenan AB, Torre D, Lachmann A, Leong AK, Wojciechowicz ML, Utti V, Jagodnik KM, Kropiwnicki E, Wang Z, Ma'ayan A. ChEA3: transcription factor enrichment analysis by orthogonal omics integration. Nucleic Acids Res. 2019 Jul 2;47(W1):W212-W224. doi: 10.1093/nar/gkz446.

