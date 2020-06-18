# KinomeScan-Appyter
Appyter to visualize KinomeScan data and integrate it with KEA3

This appyter creates bar charts for the visualization of KINOMEScan assay data from the Harvard Medical School Library of Integrated Network-based Cellular Signatures (HMS LINCS) database.


The user can input a small molecule and/or a kinase. 


If a small molecule is inputted, the appyter will return the list of kinases it binds to. The appyter will generate a bar chart displaying the kinases, sorted by their equilibrium dissociation constant, Kd (which is an indication of binding affinity). Users can hover over each bar to see the kinases associated with each equilibrium dissociation constant. 


Similarly, if a kinase is inputted, the appyter will return the list of small molecules that bind to it. The appyter will also generate an interactive bar chart displaying the small molecules sorted by their equilibrium dissociation constant. 
