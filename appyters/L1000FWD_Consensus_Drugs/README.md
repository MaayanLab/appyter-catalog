# L1000FWD Consensus Drugs

This Appyter takes up and down gene sets and returns the consensus mimicker and reverser drugs associated to them. It sends the up and down gene sets to [L1000FWD](https://amp.pharm.mssm.edu/l1000fwd/) to get the top 50 mimicker and reverser signatures. It then counts the number of times a drug appear in the top 50 signatures of a query gene set. The appyter only keeps drugs that appear in a certain percentage of query gene sets defined by **drug percentage**. The p-value for each drug is then computed using a pre-computed empirical distribution using [CREEDS](https://amp.pharm.mssm.edu/creeds/#downloads) gene sets as background signatures. The drugs are then visualized using [clustergrammer](https://amp.pharm.mssm.edu/clustergrammer/). A publication ready static heatmap is also provided for the users.

## References
[1] Wang Z, Lachmann A, Keenan AB, Ma'ayan A (2018) L1000FWD: fireworks visualization of drug-induced transcriptomic signatures. Bioinformatics doi: 10.1093/bioinformatics/bty060

[2] Wang, Z., Monteiro, C. D., Jagodnik, K. M., Fernandez, N. F., Gundersen, G. W., ... & Ma'ayan, A. (2016) Extraction and Analysis of Signatures from the Gene Expression Omnibus by the Crowd. Nature Communications doi: 10.1038/ncomms12846

[3] Fernandez, N. F. et al. Clustergrammer, a web-based heatmap visualization and analysis tool for high-dimensional biological data. Sci. Data 4:170151 doi: 10.1038/sdata.2017.151 (2017).