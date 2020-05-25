

#### upload druggable genome data

gcloud compute scp /Users/catherinestorm/Documents/UCL_PhD/PhD_Project/druggable_genome/snps_5kb_window/masterscripts/mr_druggable_genome_masterscript_test/druggable_genome_new.txt catherinestorm@pdtreatment://home/catherinestorm/





#### download the eqtl exposure data

# psychencode

mkdir eqtl_data_psychencode

wget http://resource.psychencode.org/Datasets/Derived/QTLs/DER-08a_hg19_eQTL.significant.txt -P eqtl_data_psychencode

wget http://resource.psychencode.org/Datasets/Derived/QTLs/SNP_Information_Table_with_Alleles.txt -P eqtl_data_psychencode


# eqtlgen

mkdir eqtl_data_eqtlgen

wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz -P eqtl_data_eqtlgen

wget 
https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz -P eqtl_data_eqtlgen





#### upload pd gwas outcome data

mkdir outcome_data


# upload risk data – discovery phase
gcloud compute scp AllResults_updatedrsid.txt catherinestorm@pdtreatment://home/catherinestorm/outcome_data/


# upload risk data – replication phase
gcloud compute scp toMeta* catherinestorm@pdtreatment://home/catherinestorm/outcome_data/
nohup bash script_meta_analysis_meta5_without_nalls2014.txt &
gcloud compute scp METAANALYSIS1_standarderror.TBL.zip catherinestorm@pdtreatment://home/catherinestorm/outcome_data/


# upload age at onset data
gcloud compute scp Blauwendraat_IPDGC_only_AAO_GWAS_sumstats_april_2018.txt.zip catherinestorm@pdtreatment://home/catherinestorm/outcome_data/


## download/upload the raw progression data
gcloud compute scp progression.tar.gz catherinestorm@pdtreatment://home/catherinestorm/outcome_data/

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


while read OUTCOME; do
    gcloud compute scp ${OUTCOME}.txt.gz catherinestorm@pdtreatment://home/catherinestorm/outcome_data/
done < progression_outcomes.txt



