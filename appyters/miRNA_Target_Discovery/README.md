# miRNA Target Discovery via CLIP-SEQ

This Appyter gets miRNA-gene targets from the output of dCLIP (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r11). The gene set enrichment is performed to these gene targets using Enrichr (https://amp.pharm.mssm.edu/Enrichr).


### dCLIP workflow
1. Adapter trimming and quality filtering of CLIP-SEQ reads.
2. Sequence alignment using programs like STAR or HISAT2
3. run dCLIP on the sam files see [manual here](https://qbrc.swmed.edu/download/README3.txt)