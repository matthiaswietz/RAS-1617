## "The Polar Night Shift: Annual Dynamics and Drivers of Microbial Community Structure in the Arctic Ocean"

This repo describes the amplicon-seq workflow in Wietz et al. (https://www.biorxiv.org/content/10.1101/2021.04.08.436999v2.full). Microbial samples were derived from autonomous Remote Access Samplers between July 2016 to July 2017 in different locations of the Fram Strait (Arctic Ocean). Here, we study biodiversity patterns of bacteria, archaea and microbial eukaryotes in ice-covered and ice-free regions, providing a high-resolution portrait of marine microbial ecology over polar day and night. The project is part of the long-term ecological observatory FRAM of the Alfred Wegener Institute (https://www.awi.de/en/expedition/observatories/ocean-fram.html). 

## Contents
- Rmarkdowns describing primer clipping (using cutadapt) and ASV classification (using DADA2) for [bacteria](./cutadapt_dada/cutadapt_dada_bac.Rmd) and [eukaryotes](./cutadapt_dada/cutadapt_dada_euk.Rmd)  
- [Resulting ASV and taxonomy tables](./cutadapt_dada)   
- [Rscripts and metadata](./analysisCode) to reproduce analyses and figures
