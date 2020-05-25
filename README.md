# mr_druggable_genome_pd

## Introduction
This is a collection of scripts used to complete a Mendelian randomization analysis for the druggable genome in Parkinson's disease. The method uses expression quantitative trait loci (eQTL) from blood and brain tissue for druggable genes to predict the efficacy of using these medications in Parkinson's disease. Full methods can be found [here]().
* This study used expression quantitative trait loci (eQTL) provided by the [eQTLGen](https://eqtlgen.org) and [PsychENCODE](http://resource.psychencode.org) consortia.
* This study used outcome data for PD risk (discovery and replication phase), age at onset, and progression provided by the [International Parkinson's disease Genomics consortium](http://pdgenetics.org/resources)

This code can be used for a Mendelian randomization analysis of the druggable genome using any QTL data and any disease outcome. Note!
* You will need to provide the druggable genome file. The publicly available version of the druggable genome provided by [Finan at al.](https://pubmed.ncbi.nlm.nih.gov/28356508/) can be used instead of `druggable_genome_new.txt`. ?rewrite code so it works for this
* You will need to provide the GWAS summary statistics for your outcome of interest and create `read_exposure_data_${YourExposureData}.R` and `read_exposure_data_${YourOutcomeData}.R` scripts to suit these. Use `read_exposure_data_generic.R` and `read_exposure_data_generic.R` as a template.




## Pipeline Overview



1. Install required tools

```bash
bash ./mr_druggable_genome_pd/shell/installing_tools.sh
```


2. Download required data

```bash
bash ./mr_druggable_genome_pd/shell/eqtl_data_download.sh
```


3. Select exposure data and outcome data

Edit as appropriate
```bash
echo "psychencode
eqtlgen" > exposure_data.txt
```


For discovery analyses:
```bash
echo "nalls2014
cont_HY
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
surv_DYSKINESIAS" > outcomes.txt
```


For replication analyses, make sure the outcome should begin with "replication_":
```bash
echo "replication_risk
replication_aao" > outcomes.txt
```


4. Prepare the data for the Mendelian randomization analysis. Note: if you are using your own GWAS data you will need to edit these files and create appropriate `read_exposure_data` and `read_outcome_data` files
```bash
bash ./mr_druggable_genome_pd/shell/data_prep.sh
```

5. Generate scripts that can be run in parallel
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        export DISCOVERY_OUTCOME="nalls2014" #type in
        mkdir ${EXPOSURE_DATA}_${OUTCOME}
        
        bash ./mr_druggable_genome_pd/shell/generate_parallel_scripts.sh
    done < outcomes.txt
done < exposure_data.txt
```

6. Run Mendelian randomization for the druggable genome using all the exposure data and outcome data you have specified.
```bash
echo "exposure_data,exposure,outcome" > exposures_to_remove.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./mr_druggable_genome_pd/shell/run_liberal_scripts_all_nohup.sh &> ./mr_druggable_genome_pd/shell/nohup_run_liberal_scripts_all.log &
    done < outcomes.txt
done < exposure_data.txt
```

7. Some genes cause errors during clumping or MR analysis methods using a linkage disequilibrium LD matrix. The below script will put any genes that need to be removed in `exposures_to_remove.txt` and rerun any scripts that encountered an error. Note: you may need to rerun this step a few times until all genes that cause an error are removed.
```bash
nohup bash ./mr_druggable_genome_pd/shell/run_liberal_scripts_failed_nohup.sh &> ./mr_druggable_genome_pd/shell/nohup_run_liberal_scripts_failed.log &
```

8. For each exposure-data-outcome combination, put all the results into one results file.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        cd ./mr_druggable_genome_pd/${EXPOSURE_DATA}_${OUTCOME}/results
        nohup Rscript ../R/combine_results_liberal_r2_0.2.R &> ../${EXPOSURE_DATA}_${OUTCOME}/nohup_combine_results_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.log &
        cd ..
    done < outcomes.txt
done < exposure_data.txt
```

9. As a quality control step, rerun the analysis for all genes reaching significance using a clumping threshold of 0.001.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./mr_druggable_genome_pd/${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.sh &> ./mr_druggable_genome_pd/${EXPOSURE_DATA}_${OUTCOME}/nohup_script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.log &
    done < outcomes.txt
done < exposure_data.txt
```


## Steps that are specific to this analysis
1. Create final results files and generate a report of the number of genes tested vs reaching significance for each exposure-data-outcome combination.
```bash

mkdir full_results
Rscript ./mr_druggable_genome_pd/R/final_results_report.R &> ./mr_druggable_genome_pd/nohup_final_results_report.log &

wait

cat nohup_final_results_report.log | sed "s/\[1\] //g" | sed "s/[\"]//g" > full_results/final_results_report.txt

```

2. Display the results in a forest plot


## Citation
If you use the code, please cite:
add citation
note to see methods
