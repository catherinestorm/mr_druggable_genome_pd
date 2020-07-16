
library("dplyr")
library("readr")
library("stringr")

# read in gwas data
nalls2014 <- read_tsv("outcome_data/AllResults_updatedrsid.txt")

nalls2014$chr_pos <- str_c("chr",nalls2014$Chr, ":",nalls2014$Bp)

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
nalls2014_1 <- left_join(nalls2014, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = "chr_pos")
nalls2014_2 <- left_join(nalls2014_1, psychencode[,c("chr_pos", "rsid_psychencode")], by = "chr_pos")
nalls2014_3 <- left_join(nalls2014_2, pqtl[,c("chr_pos", "rsid_replication_pqtl")], by = "chr_pos")


write.table(nalls2014_3,"outcome_data/nalls2014_discovery_risk.txt", sep = "\t", row.names = F, quote = F)

print("mission complete")
