# Enrichment Analysis Visualization Appyter

This appyter creates a variety of visualizations for enrichment analysis data for one selected Enrichr library. The scatter plot and hexagonal canvas may not be available for every library.

Here are the visualizations created: 
 
* **Scatter Plot** 

The scatterplot is organized so that simliar gene sets are clustered together. Enriched gene sets will appear blue instead of gray; the darker the blue the smaller the p-value. The name of the gene set and the associated p-value are displayed when a point is hovered over.

Plots can be downloaded as an svg using the save function on the toolbar next to the plot.

* **Bar Chart**

The bar chart contains the top 5 enriched terms and their corresponding p-values for the chosen library. Colored bars correspond to terms with significant p-values (<0.05). An asterisk (*) next to a p-value indicates the term also has a significant adjusted p-value (<0.05).

* **Hexagonal Canvas**

The hexagonal canvas is made up of hexagons, each of which represents one gene set from the library selected. The hexagons are colored based on the Jaccard similarity index between the inputted gene list and the gene set the hexagon represents (the brighter, the more similar). The hexagons representing the most similar gene sets are grouped together. The name of the gene set and the associated similarity index are displayed when a hexagon is hovered over.

* **Manhattan Plot**

A manhattan plot is created displaying all the gene sets from the select library and their p-value. The x-axis of the plot is made up of gene sets from the library. The y-axis of the plot has the -log(p value) for each gene set. The name of the gene set and the associated p-value are displayed when a point is hovered over. You can also zoom, pan, and save the plot as an svg using the toolbar on the right.

Additionally, the names and p-values of significant terms in the library are provided in a table output below the plot which can be downloaded. 

There is a link at the bottom to the full analysis results on the Enrichr website.