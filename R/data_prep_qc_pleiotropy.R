

library("dplyr")
library("readr")


dat <- read_tsv("full_results/dat_steiger_liberal_all.txt", col_types=cols(.default="c"))

sign <- read_tsv("full_results/significant_genes_results_all_outcomes.txt", col_types=cols(.default="c"))

dat_keep <- subset(dat, dat$tissue == "eqtlgen" & !(dat$outcome %in% "nalls2019"))
sign_keep <- subset(sign, sign$tissue_qtl_type == "blood_eqtl")

dat_keep1 <- subset(dat_keep, dat_keep$exposure %in% sign_keep$exposure)

# load eqtlgen data

eqtlgen <- read_tsv(gzfile("eqtl_data_eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"), col_types=cols(.default="c"))

eqtlgen$new_gene_id <- eqtlgen$Gene


# alt allele == effect allele
alleles <- read_tsv(gzfile("eqtl_data_eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"), col_types=cols(.default="c"))


# keep the genes that have at least 1 snp overlapping with a sign. gene in our data

eqtlgen_genes_with_overlapping_snps <- subset(eqtlgen, eqtlgen$SNP %in% dat_keep1$SNP)

eqtlgen_keep <- subset(eqtlgen, eqtlgen$Gene %in% eqtlgen_genes_with_overlapping_snps$Gene)

# add allele data

alleles_keep <- subset(alleles, alleles$SNP %in% eqtlgen_keep$SNP)

# AlleleB_all ==  Allele frequency of Allele B
full <- left_join(eqtlgen_keep, alleles[, c("SNP", "AlleleB", "AlleleB_all")], by = "SNP")



#switch allele frequencies where appropriate
mismatch <- which(full$Assessed != full$AlleleB)
mismatch2 <- which(full$Other == full$AlleleB)
sum(mismatch - mismatch2) # should be 0

full$eaf <- full$AlleleB_all
full$eaf[mismatch] <- 1 - as.numeric(full$eaf[mismatch])



# calculate beta and standard error

full$beta <- as.numeric(full$Zscore) / sqrt(2 * as.numeric(full$eaf) *
                                                      (1- as.numeric(full$eaf)) *
                                                      (as.numeric(full$NrSamples) + as.numeric(full$Zscore)^2))

full$se = 1 / sqrt(2 * as.numeric(full$eaf) *
                         (1- as.numeric(full$eaf)) *
                         (as.numeric(full$NrSamples) + as.numeric(full$Zscore)^2))


length(unique(full$GeneSymbol))
# add gene names

#full_with_names <- left_join(full, genes_id[, c("exposure", "gene.exposure")], by = c("new_gene_id" = "exposure"))

write.table(full, "eqtl_data_eqtlgen/eqtlgen_exposure_dat_qc.txt", sep = "\t", row.names = F, quote = F)


print("mission complete")










# psychencode

dat_keep <- subset(dat, dat$tissue == "psychencode" & !(dat$outcome %in "nalls2019"))
sign_keep <- subset(sign, sign$tissue_qtl_type == "brain_eqtl")

dat_keep1 <- subset(dat_keep, dat_keep$exposure %in% sign_keep$exposure)



# load psychencode data

psychencode <- read_tsv("eqtl_data_psychencode/DER-08a_hg19_eQTL.significant.txt", col_types=cols(.default="c"))

psychencode$new_gene_id <- gsub("\\..*", "", psychencode$gene_id)


# alt allele == effect allele
alleles <- read_tsv("eqtl_data_psychencode/SNP_Information_Table_with_Alleles.txt", col_types=cols(.default="c"))

#alleles_keep <- subset(alleles, alleles$PEC_id %in% psychencode$SNP_id)

full0 <- left_join(psychencode, alleles, by = c("SNP_id" = "PEC_id"))

# keep the genes that have at least 1 snp overlapping with a sign. gene in our data
genes_with_overlapping_snps <- subset(full0, full0$Rsid %in% dat_keep1$SNP)

full <- subset(full0, full0$new_gene_id %in% genes_with_overlapping_snps$new_gene_id)


# calculate standard error from beta and pvalue

full$se=abs(as.numeric(full$regression_slope)/qnorm(as.numeric(full$nominal_pval)/2))


# add gene names

#full_with_names <- left_join(full, genes_id[, c("exposure", "gene.exposure")], by = c("new_gene_id" = "exposure"))

length(unique(full$new_gene_id))

write.table(full, "eqtl_data_psychencode/psychencode_exposure_dat_qc.txt", sep = "\t", row.names = F, quote = F)
