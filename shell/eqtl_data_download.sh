#!/bin/bash

# Download publicly available eQTL data



### psychencode

mkdir eqtl_data_psychencode

wget http://resource.psychencode.org/Datasets/Derived/QTLs/DER-08a_hg19_eQTL.significant.txt -P eqtl_data_psychencode

wget http://resource.psychencode.org/Datasets/Derived/QTLs/SNP_Information_Table_with_Alleles.txt -P eqtl_data_psychencode



### eqtlgen

mkdir eqtl_data_eqtlgen

wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz -P eqtl_data_eqtlgen

wget 
https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz -P eqtl_data_eqtlgen

