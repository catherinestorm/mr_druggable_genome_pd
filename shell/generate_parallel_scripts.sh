#!/bin/bash

########## GENERATE SCRIPTS THAT CAN BE RUN IN PARALLEL ##########

            mkdir ${EXPOSURE_DATA}_${OUTCOME}/data
            mkdir ${EXPOSURE_DATA}_${OUTCOME}/results
            mkdir ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts


            # generate generic script for MR – clumping at r2 < 0.2
            cat ./mr_druggable_genome_pd/R/read_exposure_data_${EXPOSURE_DATA}.R > ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R

            cat ./mr_druggable_genome_pd/R/read_outcome_data_${OUTCOME}.R >> ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R

            cat ./mr_druggable_genome_pd/R/mr_analysis_liberal_r2_0.2.R >> ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R



            # generate generic script for MR with stricter clump (r2 < 0.001)
            cat ./mr_druggable_genome_pd/R/read_exposure_data_${EXPOSURE_DATA}.R > ${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.R

            cat ./mr_druggable_genome_pd/R/read_outcome_data_${OUTCOME}.R >> ${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.R

            cat ./mr_druggable_genome_pd/R/mr_analysis_conservative_r2_0.001.R >> ${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.R

    echo "export START=1
                export END='no_end'

                export EXPOSURE_DATA=${EXPOSURE_DATA}
                export OUTCOME=${OUTCOME}

                export DISCOVERY_OUTCOME=${DISCOVERY_OUTCOME}

                Rscript ${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.R" > ${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.sh






if [[ ${OUTCOME} == replication_* ]] || [[ ${EXPOSURE_DATA} == replication_* ]]; then

                echo "export START=1
                export END='no_end'

                export EXPOSURE_DATA=${EXPOSURE_DATA}
                export OUTCOME=${OUTCOME}

                export DISCOVERY_OUTCOME=${DISCOVERY_OUTCOME}

                Rscript ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R" > ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.sh


else
            if [[ $EXPOSURE_DATA == eqtlgen* ]]; then
            for ((k=1; k <= 2777;k+=50))
                do

                echo "export START=${k}
                export END=$((${k}+50))

                export EXPOSURE_DATA=${EXPOSURE_DATA}
                export OUTCOME=${OUTCOME}

                Rscript ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R" > ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh

                done
            fi

            if [[ $EXPOSURE_DATA == psychencode* ]]; then
            for ((k=1; k <= 2448;k+=50))
                do

                echo "export START=${k}
                export END=$((${k}+50))

                export EXPOSURE_DATA=${EXPOSURE_DATA}
                export OUTCOME=${OUTCOME}

                Rscript ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R" > ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh

                done
            fi

            if [[ $EXPOSURE_DATA == metabrain_bg* ]]; then
            for ((k=1; k <= 107;k+=50))
                do

                echo "export START=${k}
                export END=$((${k}+50))

                export EXPOSURE_DATA=${EXPOSURE_DATA}
                export OUTCOME=${OUTCOME}

                Rscript ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R" > ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh

                done
            fi

            if [[ $EXPOSURE_DATA == metabrain_cortex* ]]; then
            for ((k=1; k <= 1257;k+=50))
                do

                echo "export START=${k}
                export END=$((${k}+50))

                export EXPOSURE_DATA=${EXPOSURE_DATA}
                export OUTCOME=${OUTCOME}

                Rscript ${EXPOSURE_DATA}_${OUTCOME}/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.R" > ${EXPOSURE_DATA}_${OUTCOME}/${EXPOSURE_DATA}_parallel_scripts/script_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}_${k}.sh

                done
            fi
fi
