# Manhattan-Plot-Appyter

This appyter creates manhattan plots visualizing Enrichr data. The x-axis of the 
plot is sectioned by Enrichr libraries and sub-sectioned by gene sets from within
the libraries. The y-axis of the plot has the -log(p value) for each gene set.

The appyter will create two plots: one static, one dynamic. The dynamic plot should
open in a new window. The static plot will output in line with the notebook.

There is an option for a significance line to be displayed. There are also options for different 
labeling and legend styles of the plot and different color schemes.

For the static plot, and points that are above the significance line are labeled and displayed
 in a legend to the side of the plot. Static plots can be downloaded as a png, pdf, or svg.

Dynamic plots have a "hover" feature where if you mouseover any point, it will display
the gene set and p-value it represents. You can also pan, box zoom, reset the plot view,
and save the plot as an svg.
