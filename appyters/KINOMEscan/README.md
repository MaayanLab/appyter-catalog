# KinomeScan-Appyter
Appyter to visualize KinomeScan data and integrate it with KEA3

This appyter creates bar charts for the visualization of KINOMEScan assay data and Target Affinity Spectrum (TAS) values data from the Harvard Medical School Library of Integrated Network-based Cellular Signatures (HMS LINCS) database.


The user can input a small molecule and/or a kinase for the KINOMEscan and/or TAS data. 


If a small molecule is inputted for KINOMEscan, the appyter will return a table with the list of kinases it binds to and the corresponding equilibrium dissociation constant, Kd (the lower the constant, the higher the binding affinity) or % control (the lower the percentage, the higher the binding affinity). 


If a small molecule is inputted for TAS, the appyter will return the list of kinases it binds to with their approximate Kd values. The appyter will generate a bar chart displaying the kinases, sorted by Kd. 


For both KINOMEscan and TAS, users can hover over each bar in the generated bar chart to see the kinases associated with each equilibrium dissociation constant or % control. 



Similarly, if a kinase is inputted for KINOMEscan or TAS, the appyter will return either a table or list of small molecules that bind to it. The appyter will also generate an interactive bar chart displaying the small molecules sorted by their equilibrium dissociation constant or % control.
