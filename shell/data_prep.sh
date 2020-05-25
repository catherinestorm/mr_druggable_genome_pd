#!/bin/bash

# Process the GWAS data for Mendelian randomization analysis


## keep only eqtls within 5 kb of the target gene and calculate betas or SEs where needed

while read EXPOSURE_DATA; do
        nohup Rscript ./mr_druggable_genome_pd/R/data_prep_${EXPOSURE_DATA}.R &> nohup_data_prep_${EXPOSURE_DATA}.log &
done < exposure_data.txt

wait


## process discovery-phase PD risk data for the discovery phase
nohup Rscript ./mr_druggable_genome_pd/R/data_prep_nalls2014.R &> nohup_data_prep_nalls2014.log &


## process replication-phase PD risk & age at onset data
nohup Rscript ./mr_druggable_genome_pd/R/data_prep_replication.R &> nohup_data_prep_replication.log &


## process progression data
nohup Rscript ./mr_druggable_genome_pd/R/data_prep_iwaki2019.R &> nohup_data_prep_iwaki2019.log &

wait


echo "cont_HY
cont_MMSE
cont_MOCA
cont_SEADL
cont_UPDRS1_scaled
cont_UPDRS2_scaled
cont_UPDRS3_scaled
cont_UPDRS4_scaled
cont_UPDRS_scaled
surv_DEMENTIA
surv_DEPR
surv_DYSKINESIAS
surv_HY3" > progression_outcomes.txt


## generate read_outcome_data scripts for progression
while read OUTCOME; do
    cat ./mr_druggable_genome_pd/R/read_outcome_data_progression.R > ./mr_druggable_genome_pd/R/read_outcome_data_${OUTCOME}.R
done < progression_outcomes.txt
