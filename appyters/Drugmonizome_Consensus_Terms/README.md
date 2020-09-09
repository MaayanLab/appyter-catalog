# Drugmonizome Consensus Terms

This Appyter takes drug sets and returns the top consensus enrichment terms associated with them. The Appyter performs drug set enrichment analysis for all drug set libraries selected by the user that are available in Drugmonizome. The top consensus terms are determined using p-values from Drugmonizome enrichment results. The terms are then visualized using downloadable stacked bar plots and static heat maps.

## File Format
The Appyter takes drug set gmt files with the following format:
```
description 1\t\tDRUG A\tDRUG B\tDRUG C
description 2\t\tDRUG C\tDRUG D\tDRUG E
description 3\t\tDRUG B\tDRUG E\tDRUG F
description 4\t\tDRUG A\tDRUG B\tDRUG F
```

Below is a sample of an input file of drug sets:

* [Sample Drug Set](https://appyters.maayanlab.cloud/storage/Drugmonizome_Consensus/example.gmt)