#!/bin/bash

# Some genes cause errors during clumping or MR analysis methods using a linkage disequilibrium LD matrix. The below script will put any genes that need to be removed in `exposures_to_remove.txt` and rerun any scripts that encountered an error. Note: you may need to rerun this step a few times until all genes that cause an error are removed

## Make a list of all failed scripts
echo -n > failed_scripts_liberal_logs.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do

        grep -L "mission_complete" ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}*.log >> failed_scripts_liberal_logs.txt

                
    done < outcomes.txt
done < exposure_data.txt

cat failed_scripts_liberal_logs.txt | sed "s/nohup_//g" | sed "s/.log//g" > failed_scripts_liberal.txt
sed "s/log//g" -i failed_scripts_liberal.txt


wait






## Identify genes that cause an error and run the failed scripts
if [[ -s failed_scripts_liberal.txt ]] ; then
    echo "Uh oh, some scripts ran into an error. Removing genes and re-running..."
    
    
    rm *trouble_exposures*

### for errors in the clumping step
#### grab the lines with the error to find the codes for exposures we want to remove                  
while read LOGFILE; do grep -B 3 'Error: lexical error' "${LOGFILE}"; done < failed_scripts_liberal_logs.txt > clumping_trouble_exposures.txt

wait

#### extract the lines that start with "clumping", because it contains the exposure code
grep Clumping clumping_trouble_exposures.txt >> clumping_trouble_exposures1.txt

wait

#### extract only the exposure code
cat clumping_trouble_exposures1.txt | sed 's/\,.*//' | sed 's/[^ ]* //' >> clumping_trouble_exposures2.txt

wait

#### extract all instances of the exposure code from the .log files
while read LOGFILE; do grep --with-filename -f clumping_trouble_exposures2.txt "${LOGFILE}"; done < failed_scripts_liberal_logs.txt > clumping_trouble_exposures3.txt

wait

#### grab the lines that start with "Harmonising", because it contains the name of the exposure (gene)
grep Harmonising clumping_trouble_exposures3.txt >> clumping_trouble_exposures4.txt

wait

#### remove file path info but retain exposure data name
cat clumping_trouble_exposures4.txt | sed 's/\/.*log:/ /' | sed 's/_/ /' >> clumping_trouble_exposures5.txt



### for errors in the LD matrix step

#### grab the lines with the error to find the codes for exposures we want to remove                  
while read LOGFILE; do grep --with-filename -B 3 'Error in' "${LOGFILE}"; done < failed_scripts_liberal_logs.txt > ld_matrix_trouble_exposures.txt

wait

#### extract the lines that start with "exposure", because it contains the exposure code
grep exposure ld_matrix_trouble_exposures.txt > ld_matrix_trouble_exposures1.txt

wait

#### remove file path info but retain exposure data name
cat ld_matrix_trouble_exposures1.txt | sed 's/\/.*log- -//' | sed 's/_/ /' >> ld_matrix_trouble_exposures2.txt

wait

#### update the file with exposures to remove
nohup Rscript ./mr_druggable_genome_pd/R/extract_exposures_causing_trouble.R &> nohup_extract_exposures_causing_trouble.log &

wait

### rerun the failed scripts
    if [[ -s failed_scripts_liberal.txt ]] ; then
        while read FAILED; do
            export LOG=$(sed "s/script_/nohup_script_/g" <<< ${FAILED})
            nohup bash ${FAILED}.sh &> ${LOG}.log &
            wait
        done < failed_scripts_liberal.txt  
    else
        echo "Error! Is failed_scripts_liberal.txt empty?"
    fi


else
    echo "Mission complete! Move on to the next step."
fi
