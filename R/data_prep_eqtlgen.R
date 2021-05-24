

library("dplyr")
library("readr")


# load eqtlgen data

eqtlgen <- read_tsv(gzfile("eqtl_data_eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"), col_types=cols(.default="c"))

eqtlgen$new_gene_id <- eqtlgen$Gene


# alt allele == effect allele
alleles <- read_tsv(gzfile("eqtl_data_eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"), col_types=cols(.default="c"))


# load druggable genome

druggable <- read_tsv("druggable_genome_new.txt", col_types=cols(.default="c"))


# keep eqtls for druggable genes

eqtlgen_druggable <- subset(eqtlgen, eqtlgen$new_gene_id %in% druggable$gene_stable_id)







### select SNPs within 5kb upstream/downstream of gene


# read in genes positions from druggable genome

genes_id0 <- distinct(druggable[, c("gene_stable_id", "gene_display_label", "chr_name", "gene_start", "gene_end")])



# keep autosome only

genes_id01 <- genes_id0[which(genes_id0$chr_name %in% 1:22),]

genes_id01[,3:5] <- sapply(genes_id01[,3:5], as.numeric)



# renames columns

names(genes_id01)<-c("exposure","gene.exposure","chromosome_name","start_position","end_position")



# keep genes that have an eqtl eqtlgen data

genes_id <- subset(genes_id01, genes_id01$exposure %in% eqtlgen$new_gene_id)




# format eqtlgen to suit loop

dat <- eqtlgen_druggable

dat[, c("SNPChr", "SNPPos")] <- sapply(dat[, c("SNPChr", "SNPPos")], as.numeric)


# create empty data frame

genes_data <- data.frame()


# loop to keep snps within 5kb of gene start/end positions # from 1382626 eQTLs for 2918 genes to 280216 eQTLs for 2791 genes

for (i in 1:length(unique(genes_id$gene.exposure))) {

  dat1 <- dat[which(dat$new_gene_id==genes_id$exposure[i] & dat$SNPChr==genes_id$chromosome_name[i] & dat$SNPPos >= (genes_id$start_position[i]-5000) & dat$SNPPos <= (genes_id$end_position[i]+5000)),]

  genes_data <- rbind(genes_data,dat1)

}


# add allele data

alleles_keep <- subset(alleles, alleles$SNP %in% genes_data$SNP)

# AlleleB_all ==  Allele frequency of Allele B
full <- left_join(genes_data, alleles[, c("SNP", "AlleleB", "AlleleB_all")], by = "SNP")



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



# add gene names

full_with_names <- left_join(full, genes_id[, c("exposure", "gene.exposure")], by = c("new_gene_id" = "exposure"))

write.table(full, "eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t", row.names = F, quote = F)


print("mission complete")
