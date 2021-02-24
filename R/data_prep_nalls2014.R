
library("dplyr")
library("readr")
library("stringr")

# read in gwas data
pd_risk_discovery <- read_tsv("outcome_data/AllResults_updatedrsid.txt")

pd_risk_discovery$chr_pos <- str_c("chr",pd_risk_discovery$Chr, ":",pd_risk_discovery$Bp)

# read in eqtl data
eqtlgen <- read.table("eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c("chr",eqtlgen$chr, ":",eqtlgen$position)


psychencode <- read.table("eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

psychencode <- distinct(psychencode[,c("Rsid", "chr", "position")])

names(psychencode) <- c("rsid_psychencode", "chr", "position")

psychencode$chr_pos <- str_c(psychencode$chr, ":",psychencode$position)





# pqtl
pqtl <- read.table("pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt", sep = "\t",colClasses = "character", header = T)

pqtl$chr_pos <- str_c("chr", pqtl$chr, ":", pqtl$pos)
names(pqtl)[1] <- "rsid_replication_pqtl"

# join to add rsids to the gwas data
pd_risk_discovery_1 <- left_join(pd_risk_discovery, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = "chr_pos")
pd_risk_discovery_2 <- left_join(pd_risk_discovery_1, psychencode[,c("chr_pos", "rsid_psychencode")], by = "chr_pos")
pd_risk_discovery_3 <- left_join(pd_risk_discovery_2, pqtl[,c("chr_pos", "rsid_replication_pqtl")], by = "chr_pos")


write.table(pd_risk_discovery_3,"outcome_data/pd_risk_discovery_discovery_risk.txt", sep = "\t", row.names = F, quote = F)

print("mission complete")
