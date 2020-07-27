# Drugmonizome ETL: STITCH

[STITCH](http://stitch.embl.de/) is a database of known and predicted interactions between chemicals and proteins. The interactions include direct (physical) and indirect (functional) associations; they stem from computational prediction, from knowledge transfer between organisms, and from interactions aggregated from other (primary) databases.

This appyter takes drug-protein associations from STITCH and outputs files that are usable for Drugmonizome. These data are processed in order to construct a binary matrix with small molecules as rows and genes as column attributes. From this matrix, drug and attribute set libraries are constructed.

Additionally, gene and attribute similarity matrices, which store the jaccard distance between any two genes or attributes, are created.

The downloadable file will have the following outputs:
* Edge list of drug-attribute associations
* Binary matrix of drug-attribute associations
* Drug set library: matches each attribute with a set of associated small molecules
* Attribute set library: matches each small molecule with a set of associated attributes
* Drug similarity matrix
* Attribute similarity matrix