# ChEA-KG-TS Appyter

Given either raw time series RNA-seq data or pre-computed DEGs at each time point, the ChEA-KG-TS Appyter visualizes on a regulatory subnetwork how the enriched TFs at one time point regulate the enriched TFs at the subsequent time point, therefore enabling users to understand how the TF landscape governing a biological process evolves over time. The Appyter also determines which modules of TFs may be targeting similar genes at each time point by visualizing how the enriched TFs cluster on a UMAP plot. The Appyter also outputs the average rank of each enriched TF within the ChEA3 gene set libraries, enabling users to see which TFs may play a more substantial role in dictating the gene expression changes at each time point.

## Launching the appyter
```bash
appyter --profile=biojupies chea_kg_ts_appyter.ipynb
```

As long as the terminal is open, the appyter should be available at https://localhost:5000/

## Developing the appyter
The appyter can be developed with any jupyter notebook compatible editor.
```bash
jupyter chea_kg_ts_appyter.ipynb
```
