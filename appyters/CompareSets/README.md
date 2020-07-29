# Compare-Sets-Appyter
This appyter visualizes the intersections between 2-6 user-inputted gene/drug sets. The user has a choice whether to upload a tsv file with the desired genes or to insert the genes into text boxes. Three default gene sets are provided. 

Visualization options include Venn diagrams, SuperVenn diagrams, and UpSet plots. UpSet plots can be more useful than Venn diagrams when trying to visualize more than three or four sets. SuperVenn diagrams can be useful when the user wants to visualize the intersections in a tabular format. Users can customize their figures by changing the colors and other features. 

The items in the set intersections are displayed in a table with links to Enrichr for further enrichment analysis.

Fisher's Exact Test is also computed to see whether the overlap of two sets is significant. After adding their own sets, the user can select their desired significance level and background, and Fisher's Exact Test will be calculated for all pairs of sets. This data is displayed in a table with a dropdown menu to select which overlap the user would like to examine. A heatmap is also generated to display the results of all Fisher's Exact Test calculations. 
