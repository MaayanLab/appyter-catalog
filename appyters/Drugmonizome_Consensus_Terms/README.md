# Drugmonizome Consensus Terms

This Appyter takes as input a collection of drug sets and returns the top consensus enrichment terms associated with them. The Appyter performs drug set enrichment analysis for all uploaded drug sets against drug set libraries that the user selects from all the libraries that are available in Drugmonizome. The top consensus terms are determined using p-values returned from the Drugmonizome enrichment results using the Fisher exact test. The enriched terms are then visualized using downloadable stacked bar plots and static heatmaps.

## File Format
The Appyter takes as input a collection of drug sets organized in the following format:
```
description 1\t\tDRUG A\tDRUG B\tDRUG C
description 2\t\tDRUG C\tDRUG D\tDRUG E
description 3\t\tDRUG B\tDRUG E\tDRUG F
description 4\t\tDRUG A\tDRUG B\tDRUG F
```

Below is a sample of an input file containing a collection of drug sets:

* [Sample Collection of Drug Sets](https://appyters.maayanlab.cloud/storage/Drugmonizome_Consensus/example.gmt)