
library("dplyr")
library("readr")
library("stringr")

# read in allele data
alleles <- read.table(gzfile("outcome_data/reference.txt.gz"), sep = ",", stringsAsFactors = F, header = T)
names(alleles) <- c("chrpos", "rsid", "chr", "start", "other_allele", "effect_allele", "maf", "func", "near_gene")

# read in QTL data

# eqtlgen
eqtlgen <- read.table("eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c(eqtlgen$chr, ":",eqtlgen$position)

## psychencode
psychencode <- read.table("eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

psychencode <- distinct(psychencode[,c("Rsid", "chr", "position")])

names(psychencode) <- c("rsid_psychencode", "chr", "position")

psychencode$chr_pos <- str_c(psychencode$chr, ":",psychencode$position)
psychencode$chr_pos <- gsub("chr", "", psychencode$chr_pos)

## pqtl
pqtl <- read.table("pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt", sep = "\t",colClasses = "character", header = T)

pqtl$chr_pos <- str_c(pqtl$chr, ":", pqtl$pos)
names(pqtl)[1] <- "rsid_replication_pqtl"


# read in  each file and populate a data frame
file_list_res_a <- list.files("outcome_data", pattern = "cont_", full.names = T)
file_list_res_b <- list.files("outcome_data", pattern = "surv_", full.names = T)
file_list_res <- c(file_list_res_a, file_list_res_b)


for (i in 1:length(file_list_res)) {
  raw_data <- read.table(gzfile(file_list_res[i]), sep = "\t", stringsAsFactors = F, header = T)

  raw_data$outcome <- file_list_res[i]
  raw_data$outcome <- gsub(".txt.gz", "", raw_data$outcome)
  raw_data$outcome <- gsub("outcome_data/", "", raw_data$outcome)

  raw_data_with_alleles <- left_join(raw_data, alleles[, c("chrpos", "rsid", "other_allele", "effect_allele", "maf")], by = c("SNP" = "chrpos"))

  raw_data_with_alleles_1 <- left_join(raw_data_with_alleles, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = c("SNP" = "chr_pos"))
  raw_data_with_alleles_1 <- left_join(raw_data_with_alleles_1, psychencode[,c("chr_pos", "rsid_psychencode")], by = c("SNP" = "chr_pos"))
  raw_data_with_alleles_1 <- left_join(raw_data_with_alleles_1, pqtl[,c("chr_pos", "rsid_replication_pqtl")], by = c("SNP" = "chr_pos"))

  write.table(raw_data_with_alleles_1, str_c("outcome_data/", paste(raw_data_with_alleles$outcome[1]), "_with_alleles.txt"), sep = "\t", row.names = F, quote = F)
}

print("mission complete")
