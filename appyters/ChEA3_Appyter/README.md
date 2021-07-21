# Example

This is an example appyter in the catalog. Please reference it when submitting your own appyter.

- `{name}/{name}.ipynb`: ChEA3 Appyter
- `{name}/README.md`: The ChEA3 Appyter (ChIP-X Enrichment Analysis 3) predicts transcription factors (TFs) associated with user-input sets of genes. Discrete query gene sets are compared to ChEA3 libraries of TF target gene sets assembled from multiple orthogonal 'omics' datasets. The Fisher's Exact Test, with a background size of 20,000, is used to compare the input gene set to the TF target gene sets in order to determine which TFs may be most closely associated with the input gene set. Users can submit a set of human or mouse gene symbols for transcription factor enrichment analysis.
- `{name}/deps.txt`: Docker image (ubuntu) packages necessary to run your appyter
- `{name}/requirements.txt`: requests, numpy, tabulate, plotly, kaleido
- `{name}/appyter.json`: Metadata about your appyter to be used for the website, be sure to make `"public": true` if you want it visible in the catalog.

![Screenshot](./static/screenshot.png)
