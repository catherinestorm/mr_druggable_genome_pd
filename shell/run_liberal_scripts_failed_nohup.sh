rm *trouble_exposures*

### for errors in the clumping step
# grab the lines with the error to find the codes for exposures we want to remove                  
while read LOGFILE; do grep -B 3 'Error: lexical error' "${LOGFILE}"; done < failed_scripts_liberal_logs.txt > clumping_trouble_exposures.txt

wait

# extract the lines that start with "clumping", because it contains the exposure code
grep Clumping clumping_trouble_exposures.txt >> clumping_trouble_exposures1.txt

wait

# extract only the exposure code
cat clumping_trouble_exposures1.txt | sed 's/\,.*//' | sed 's/[^ ]* //' >> clumping_trouble_exposures2.txt

wait

# extract all instances of the exposure code from the .log files
while read LOGFILE; do grep --with-filename -f clumping_trouble_exposures2.txt "${LOGFILE}"; done < failed_scripts_liberal_logs.txt > clumping_trouble_exposures3.txt

wait

# grab the lines that start with "Harmonising", because it contains the name of the exposure (gene)
grep Harmonising clumping_trouble_exposures3.txt >> clumping_trouble_exposures4.txt

wait

# remove file path info but retain exposure data name
cat clumping_trouble_exposures4.txt | sed 's/\/.*log:/ /' | sed 's/_/ /' >> clumping_trouble_exposures5.txt



### for errors in the LD matrix step

# grab the lines with the error to find the codes for exposures we want to remove                  
while read LOGFILE; do grep --with-filename -B 3 'Error in' "${LOGFILE}"; done < failed_scripts_liberal_logs.txt > ld_matrix_trouble_exposures.txt

wait

# extract the lines that start with "exposure", because it contains the exposure code
grep exposure ld_matrix_trouble_exposures.txt > ld_matrix_trouble_exposures1.txt

wait

# remove file path info but retain exposure data name
cat ld_matrix_trouble_exposures1.txt | sed 's/\/.*log- -//' | sed 's/_/ /' >> ld_matrix_trouble_exposures2.txt

wait

# update the file with exposures to remove
nohup Rscript extract_exposures_causing_trouble.R &> nohup_extract_exposures_causing_trouble.log &

wait

# rerun the failed scripts
if [[ -s failed_scripts_liberal.txt ]] ; then
    while read FAILED; do
        export LOG=$(sed "s/script_/nohup_script_/g" <<< ${FAILED})
        nohup bash ${FAILED}.sh &> ${LOG}.log &
        wait
    done < failed_scripts_liberal.txt  
else
    echo "Error! Is failed_scripts_liberal.txt empty?"
fi

