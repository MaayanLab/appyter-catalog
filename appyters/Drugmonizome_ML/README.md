# Drugmonizome-ML

This appyter creates a custom machine learning pipeline on top of the [Drugmonizome](https://amp.pharm.mssm.edu/drugmonizome/) drug set libraries in order to predict drug candidates for repurposing studies.

The drug set libraries associate thousands of drugs with known attributes such as associated genes, indications, side effects, and chemical features, which can be used as inputs to the classifier. The target classes are provided as a list of drugs to be labelled as positives (e.g. the hits discovered from an *in vitro* drug screen). A model is then trained on this binary classification task.

Afterwards, the trained model is used to rank the input drug list, and to predict additional drugs that could could be similar candidates to the provided drugs.

*Note:* The default list of drugs contains the combined hits from 7 *in vitro* [COVID-19 drug screens](https://amp.pharm.mssm.edu/covid19/):  
[1] Ellinger, B., and Andrea Zaliani. "Identification of inhibitors of SARS-CoV-2 in-vitro cellular toxicity in human (Caco-2) cells using a large scale drug repurposing collection." Research Square 10 (2020).  
[2] Heiser, Katie, et al. "Identification of potential treatments for COVID-19 through artificial intelligence-enabled phenomic analysis of human cells infected with SARS-CoV-2." bioRxiv (2020).  
[3] Jeon, Sangeun, et al. "Identification of antiviral drug candidates against SARS-CoV-2 from FDA-approved drugs." Antimicrobial Agents and Chemotherapy (2020).  
[4] Ko, Meehyun, et al. "Screening of FDA-approved drugs using a MERS-CoV clinical isolate from South Korea identifies potential therapeutic options for COVID-19." bioRxiv (2020).  
[5] Mirabelli, Carmen, et al. "Morphological Cell Profiling of SARS-CoV-2 Infection Identifies Drug Repurposing Candidates for COVID-19." bioRxiv (2020).  
[6] Riva, Laura, et al. "A Large-scale Drug Repositioning Survey for SARS-CoV-2 Antivirals." bioRxiv (2020).  
[7] Touret, Franck, et al. "In vitro screening of a FDA approved chemical library reveals potential inhibitors of SARS-CoV-2 replication." BioRxiv (2020).
