# mr_druggable_genome_pd

## Introduction
This is a collection of scripts used in a Mendelian randomization analysis for the druggable genome in Parkinson's disease. Full methods can be found [here](https://www.biorxiv.org/content/10.1101/2020.07.24.208975v1).

This code can be used for any QTL data and any disease outcome.
* You will need to provide the druggable genome file. The publicly available version of the druggable genome provided by [Finan at al.](https://pubmed.ncbi.nlm.nih.gov/28356508/) can be used instead of `druggable_genome_new.txt`.
* You will need to provide the GWAS summary statistics for your outcome of interest and create `read_exposure_data_${YOUR_EXPOSURE_DATA}.R` and `read_exposure_data_${YOUR_OUTCOME_DATA}.R` scripts to suit these. Use `read_exposure_data_generic.R` and `read_exposure_data_generic.R` as a template.



## Citation
If you use the code, please cite:
[Storm CS, Kia DA, Almramhi M, Bandres-Ciga S, Finan C, Hingorani AD, International Parkinson’s Disease Genomics Consortium (IPDGC), Wood NW. "Finding drug targeting mechanisms with genetic evidence for Parkinson’s disease." bioRxiv. 2020.](https://www.biorxiv.org/content/10.1101/2020.07.24.208975v1)



## Pipeline Overview



1. Install required tools

```bash
bash ./mr_druggable_genome_pd/shell/installing_tools.sh
```


2. Download/upload your QTL and disease GWAS data. eQTL data from the PsychENCODE and eQTLGen consortia can be accessed as shown below. If you are looking for a comprehensive list of pQTL data available, check out [this very helpful blog post](http://www.metabolomix.com/a-table-of-all-published-gwas-with-proteomics/).

```bash
bash ./mr_druggable_genome_pd/shell/eqtl_data_download.sh
```


3. Select exposure data and outcome data.

Specify your exposure data. Example below.
```bash
echo "psychencode
eqtlgen" > exposure_data.txt
```


For discovery analyses, specify your outcomes. Example below.
```bash
echo "pd_risk_discovery
pd_progression_cont_HY
pd_progression_cont_MMSE
pd_progression_cont_MOCA
pd_progression_cont_SEADL
pd_progression_cont_UPDRS1_scaled
pd_progression_cont_UPDRS2_scaled
pd_progression_cont_UPDRS3_scaled
pd_progression_cont_UPDRS4_scaled
pd_progression_cont_UPDRS_scaled
pd_progression_surv_DEMENTIA
pd_progression_surv_DEPR
pd_progression_surv_DYSKINESIAS" > outcomes.txt
```


For replication analyses, the outcome should begin with "replication_", and you need to specify your discovery outcome. This will make sure your pvalue threshold for significance is the non-adjusted pvalue, and that you only run the replication step for discovered genes. Examples below.
```bash
echo "replication_risk
replication_aao" > outcomes.txt

echo "replication_pqtl" > outcomes.txt

echo "pd_risk_discovery" > discovery_outcomes.txt
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

7. Remove genes that throw errors. Note: you may need to rerun this step a few times until all genes that cause an error are removed.
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

9. Rerun the analysis for all significant genes using the clumping threshold r2 = 0.001.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./mr_druggable_genome_pd/${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.sh &> ./mr_druggable_genome_pd/${EXPOSURE_DATA}_${OUTCOME}/nohup_script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.log &
    done < outcomes.txt
done < exposure_data.txt
```

10. Some data formatting steps.
```bash
mkdir full_results

bash ./mr_druggable_genome_pd/shell/final_results_report.sh

Rscript ./mr_druggable_genome_pd/R/combine_dat_steiger.R

Rscript ./mr_druggable_genome_pd/R/format_supplement.R

```


11. Check if the direction of effect is consistent between any discovery cohorts and replication cohorts.


```bash
while read OUTCOME; do
    export OUTCOME=${OUTCOME}
    export DISCOVERY_OUTCOME=${DISCOVERY_OUTCOME}

    nohup Rscript ./mr_druggable_genome_pd/R/check_direction_of_effect.R &> full_results/metric_check_direction_of_effect_${OUTCOME}_${DISCOVERY_OUTCOME}.log
done < outcomes.txt
```


13. Generate forest plots. This script is not generic for any exposure/outcome data.

```bash
mkdir figures

Rscript ./mr_druggable_genome_pd/R/make_forest_plots.R
```

14. Colocalization analysis.

```bash
mkdir coloc

Rscript ./mr_druggable_genome_pd/R/data_prep_coloc.R

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        Rscript ./mr_druggable_genome_pd/R/run_coloc.R
    done < outcomes.txt
done < exposure_data.txt

```
