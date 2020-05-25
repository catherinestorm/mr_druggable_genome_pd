

library("dplyr")
library("readr")


# load psychencode data

psychencode <- read_tsv("eqtl_data_psychencode/DER-08a_hg19_eQTL.significant.txt", col_types=cols(.default="c"))

psychencode$new_gene_id <- gsub("\\..*", "", psychencode$gene_id)


# alt allele == effect allele
alleles <- read_tsv("eqtl_data_psychencode/SNP_Information_Table_with_Alleles.txt", col_types=cols(.default="c"))



# load druggable genome

druggable <- read_tsv("druggable_genome_new.txt", col_types=cols(.default="c"))


# keep eqtls for druggable genes # 3645 genes available

psychencode_druggable <- subset(psychencode, psychencode$new_gene_id %in% druggable$gene_stable_id)







### select SNPs within 5kb upstream/downstream of gene


# read in genes positions from druggable genome

genes_id0 <- distinct(druggable[, c("gene_stable_id", "gene_display_label", "chr_name", "gene_start", "gene_end")])



# keep autosome only

genes_id01 <- genes_id0[which(genes_id0$chr_name %in% 1:22),]

genes_id01[,3:5] <- sapply(genes_id01[,3:5], as.numeric)



# renames columns

names(genes_id01)<-c("exposure","gene.exposure","chromosome_name","start_position","end_position") 



# keep genes that have an eqtl psychencode data

genes_id <- subset(genes_id01, genes_id01$exposure %in% psychencode$new_gene_id)




# format psychencode to suit loop

dat <- psychencode_druggable

dat$SNP_chr <- gsub("chr", "",dat$SNP_chr)

dat[, c("SNP_chr", "SNP_start")] <- sapply(dat[, c("SNP_chr", "SNP_start")], as.numeric)



# create empty data frame

genes_data <- data.frame()


# loop to keep snps within 5kb of gene start/end positions # from 343084 eQTLs for 3645 genes to 80069 eQTLs for 2449 genes

for (i in 1:length(unique(genes_id$gene.exposure))) {
  
  dat1 <- dat[which(dat$new_gene_id==genes_id$exposure[i] & dat$SNP_chr==genes_id$chromosome_name[i] & dat$SNP_start >= (genes_id$start_position[i]-5000) & dat$SNP_start <= (genes_id$end_position[i]+5000)),]

  genes_data <- rbind(genes_data,dat1)
  
}


# add allele data

alleles_keep <- subset(alleles, alleles$PEC_id %in% genes_data$SNP_id)

full <- left_join(genes_data, alleles, by = c("SNP_id" = "PEC_id"))


# calculate standard error from beta and pvalue

full$se=abs(as.numeric(full$regression_slope)/qnorm(as.numeric(full$nominal_pval)/2))


# add gene names

full_with_names <- left_join(full, genes_id[, c("exposure", "gene.exposure")], by = c("new_gene_id" = "exposure"))


# sanity check  # should be TRUE for all

all((as.numeric(full_with_names$SNP_distance_to_TSS)-5000) < full_with_names$SNP_start)



write.table(full_with_names, "eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt", sep = "\t", row.names = F, quote = F)


print("mission complete")

