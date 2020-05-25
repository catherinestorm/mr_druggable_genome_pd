# mr_druggable_genome_pd

## Introduction
This is a collection of scripts used to complete a Mendelian randomization analysis for the druggable genome in Parkinson's disease. The method uses expression quantitative trait loci (eQTL) from blood and brain tissue for druggable genes to predict the efficacy of using these medications in Parkinson's disease.


This code can be used for a Mendelian randomization analysis of the druggable genome using any QTL data and any disease outcome.

For full methods, please see...

## Pipeline Overview

1. Install required tools

`bash ./mr_druggable_genome_pd/shell/installing_tools.sh`


2. Download required data

Notes
* You will need to provide the druggable genome file. The publicly available version of the druggable genome provided by [Finan at al.](https://pubmed.ncbi.nlm.nih.gov/28356508/) can be used instead of `druggable_genome_new.txt`. ?rewrite code so it works for this
* You will need to provide the GWAS summary statistics for your outcome.
* This study used outcome data for PD risk (discovery and replication phase), age at onset, and progression provided by the [International Parkinson's disease Genomics consortium](http://pdgenetics.org/resources)

`bash ./mr_druggable_genome_pd/shell/eqtl_data_download.sh`





## Citation
If you use the code, please cite:
add citation
note to see methods
