# Dataset Pair Attribute Similarity

This appyter compares the similarity between two gene set libraries.

The appyter creates two attribute similarity matrices from these libraries using two set similarity measures:

* The [Jaccard Index](https://en.wikipedia.org/wiki/Jaccard_index) measures the degree of similarity between each gene set by comparing the sets' intersection over their union.
* The Fisher's Exact Test is a proportions test that measures the probability that two sets have a statistically significant overlap.

The appyter then creates: a gene set length histogram, a clustered heatmap of the Jaccard distance matrix, and a clustered heatmap of the Fisher's exact test distance matrix.

Additionally, two tables are created showing the top 20 attribute pairs from the Jaccard indices and Fisher exact test p-values distance matruces. These attribute pairs have the highest degree of overlapping sets.

Finally, two download links are provided, one for each attribute similarity matrix.
