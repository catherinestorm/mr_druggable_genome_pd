#liftover_metabrain.sh

export STUDY4LIFTOVER="eqtl_data_metabrain/metabrain_4liftover_bed"

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# Convert sum states file to bed bile for crossmap
Rscript ./mr_druggable_genome_pd/R/sumstats2bed_4crossmap.R

# Now run CrossMap
.local/bin/CrossMap.py bed hg38ToHg19.over.chain.gz ${STUDY4LIFTOVER}.txt > ${STUDY4LIFTOVER}_crossmap_hg38ToHg19.txt

# Reformat the output to be readable by R

sed 's/Unmap/fail\tfail\tfail\tfail/g' ${STUDY4LIFTOVER}_crossmap_hg38ToHg19.txt > ${STUDY4LIFTOVER}_crossmap_hg38ToHg19_readable.txt

sed -e 's/chr//g' -e 's/X/23/g' -e 's/Y/24/g' -e 's/M/26/g' ${STUDY4LIFTOVER}_crossmap_hg38ToHg19_readable.txt > ${STUDY4LIFTOVER}_crossmap_hg38ToHg19_readable_chrplink.txt

grep -v "fail" ${STUDY4LIFTOVER}_crossmap_hg38ToHg19_readable_chrplink.txt > ${STUDY4LIFTOVER}_crossmap_hg38ToHg19_readable_chrplink_mapped.txt


# Then view results in R and generate new bim file
Rscript ./mr_druggable_genome_pd/R/update_metabrain_pos_crossmap.R
