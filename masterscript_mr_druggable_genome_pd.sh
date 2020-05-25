########### data download/upload ###########

bash data_download_upload.sh &> nohup_data_download_upload.log &

########### data prep & generate scripts that can be run in parallel ###########

# edit as appropriate
echo "psychencode
eqtlgen" > exposure_data.txt


# uncomment/edit for discovery
#echo "nalls2014" > outcomes.txt
#cat progression_outcomes.txt >> outcomes.txt


#  uncomment for replication
echo "nalls2014
replication_risk
replication_aao" > outcomes.txt



### data prep

# keep only eqtls within 5 kb of the target gene and calculate betas or SEs, where appropriate

while read EXPOSURE_DATA; do
        nohup Rscript data_prep_${EXPOSURE_DATA}.R &> nohup_data_prep_${EXPOSURE_DATA}.log &
done < exposure_data.txt

wait


# process discrovery phase risk data
nohup Rscript data_prep_nalls2014.R &> nohup_data_prep_nalls2014.log &


# process replication phase risk & age at onset data
nohup Rscript data_prep_replication.R &> nohup_data_prep_replication.log &


# process progression data
nohup Rscript data_prep_iwaki2019.R &> nohup_data_prep_iwaki2019.log &

wait

# generate read_outcome_data scripts for progression
while read OUTCOME; do
    cat read_outcome_data_progression.R > read_outcome_data_${OUTCOME}.R
done < progression_outcomes.txt







### generate scripts that can be run in parallel

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        export DISCOVERY_OUTCOME="nalls2014" #type in
        mkdir ${EXPOSURE_DATA}_${OUTCOME}
        
        nohup bash generate_parallel_scripts.sh &> ${EXPOSURE_DATA}_${OUTCOME}/nohup_generate_parallel_scripts_${EXPOSURE_DATA}_${OUTCOME}.log &
    done < outcomes.txt
done < exposure_data.txt


wait








########### run liberal-clumping scripts ###########

## Step 1. run all the liberal-clumping scripts

#echo "exposure_data,exposure,outcome" > exposures_to_remove.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash run_liberal_scripts_all_nohup.sh &> nohup_run_liberal_scripts_all.log &
    done < outcomes4.txt
done < exposure_data1.txt




## Step 2. to run all the FAILED liberal-clumping scripts
# NOTE: you may need to run step 2 several times

# Step 2a. make a list of all failed scripts
echo -n > failed_scripts_liberal_logs.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do

        grep -L "mission_complete" ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}*.log >> failed_scripts_liberal_logs.txt

                
    done < outcomes.txt
done < exposure_data.txt


cat failed_scripts_liberal_logs.txt | sed "s/nohup_//g" | sed "s/.log//g" > failed_scripts_liberal.txt
sed "s/log//g" -i failed_scripts_liberal.txt

wc -l failed_scripts_liberal_logs.txt
head failed_scripts_liberal_logs.txt


wait




# Step 2c. rerun the failed scripts
if [[ -s failed_scripts_liberal.txt ]] ; then
    echo "Uh oh, some scripts ran into an error. Removing genes and re-running..."
    nohup bash run_liberal_scripts_failed_nohup.sh &> nohup_run_liberal_scripts_failed.log &
else
    echo "Mission complete! Move on to the next step."
fi





## Step 3. put all the results in one large file per exposure-data-outcome combination

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        cd ${EXPOSURE_DATA}_${OUTCOME}/results
        nohup Rscript ~/combine_results_liberal_r2_0.2.R &> ~/${EXPOSURE_DATA}_${OUTCOME}/nohup_combine_results_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.log &
        cd
    done < outcomes.txt
done < exposure_data.txt







########### run stricter-clumping scripts for quality control ###########

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.sh &> ${EXPOSURE_DATA}_${OUTCOME}/nohup_script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.log &
    done < outcomes.txt
done < exposure_data.txt




while read EXPOSURE_DATA; do
    while read OUTCOME; do
        tail ${EXPOSURE_DATA}_${OUTCOME}/nohup_script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.log 
    done < outcomes.txt
done < exposure_data.txt






########### full results report ###########

mkdir full_results
Rscript final_results_report.R &> nohup_final_results_report.log &

wait

cat nohup_final_results_report.log | sed "s/\[1\] //g" | sed "s/[\"]//g" > full_results/final_results_report.txt
less full_results/final_results_report.txt

