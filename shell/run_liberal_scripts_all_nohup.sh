#!/bin/bash

# RUNNING ALL LIBERAL SCRIPTS

if [[ ${OUTCOME} == replication_* ]] || [[ ${EXPOSURE_DATA} == replication_* ]]; then

                nohup bash ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.sh &> ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.log &


else

            if [[ $EXPOSURE_DATA == eqtlgen* ]]; then
            for ((k=1; k <= 2786;k+=50))
                do

                nohup bash ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh &> ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.log &
                wait
                done
            fi

            if [[ $EXPOSURE_DATA == psychencode* ]]; then
            for ((k=1; k <= 2448;k+=50))
                do

                nohup bash ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh &> ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.log &
                wait
                done
            fi

            if [[ $EXPOSURE_DATA == metabrain_bg* ]]; then
            for ((k=1; k <= 107;k+=50))
                do

                nohup bash ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh &> ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.log &

                done
            fi

            if [[ $EXPOSURE_DATA == metabrain_cortex* ]]; then
            for ((k=1; k <= 1257;k+=50))
                do

                nohup bash ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh &> ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/nohup_script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.log &

                wait

                done
            fi

fi
