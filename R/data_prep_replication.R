library("readr")
library("dplyr")
library("stringr")





# read in eqtl data – will get rs IDs for the SNPs from here
eqtlgen <- read.table("eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c("chr",eqtlgen$chr, ":",eqtlgen$position)


psychencode <- read.table("eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

psychencode <- distinct(psychencode[,c("Rsid", "chr", "position")])

names(psychencode) <- c("rsid_psychencode", "chr", "position")

psychencode$chr_pos <- str_c(psychencode$chr, ":",psychencode$position)



# read in replication risk data
replication_data_risk <- read_tsv("outcome_data/METAANALYSIS_samplesize_1.tbl")

replication_data_risk$beta <- replication_data_risk$Zscore / sqrt(2 * replication_data_risk$Freq1 * (1- replication_data_risk$Freq1) * (replication_data_risk$Weight + replication_data_risk$Zscore^2))
                                                  
replication_data_risk$se <- 1 / sqrt(2 * replication_data_risk$Freq1 * (1- replication_data_risk$Freq1) * (replication_data_risk$Weight + (replication_data_risk$Zscore)^2))


# join to add rsids to the gwas data
replication_data_risk <- left_join(replication_data_risk, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = c("MarkerName" = "chr_pos"))
replication_data_risk <- left_join(replication_data_risk, psychencode[,c("chr_pos", "rsid_psychencode")], by = c("MarkerName" = "chr_pos"))

write.table(replication_data_risk, "outcome_data/pd_replication_replication_risk.txt", sep = "\t", row.names = F, quote = F)



# read in replication age at onset data
replication_data_aao <- read_tsv(unzip("outcome_data/Blauwendraat_IPDGC_only_AAO_GWAS_sumstats_april_2018.txt.zip"))

replication_data_aao <- left_join(replication_data_aao, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = c("MarkerName" = "chr_pos"))
replication_data_aao <- left_join(replication_data_aao, psychencode[,c("chr_pos", "rsid_psychencode")], by = c("MarkerName" = "chr_pos"))


write.table(replication_data_aao, "outcome_data/pd_age_at_onset_replication_aao.txt", sep = "\t", row.names = F, quote = F)



print("mission complete")