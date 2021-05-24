#update_metabrain_pos_crossmap.R


library("dplyr")
library("readr")

study <- Sys.getenv("STUDY4LIFTOVER")

hg19_bed <- read.table(paste(study, "_crossmap_hg38ToHg19_readable_chrplink_mapped.txt", sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)

hg19_bed$V4 <- NULL
hg19_bed$V2 <- NULL
hg19_bed$V6 <- NULL
hg19_bed$key <- paste(hg19_bed$V1, hg19_bed$V3, sep=":")
names(hg19_bed)[1:4] <- c("original_chr", "original_pos", "new_chromosome", "new_position")


# rename contigs
#test <- hg19_bed[grepl("_g",hg19_bed$new_chromosome),]
#hg19_bed$new_chromosome <- gsub("\\_g.*","",hg19_bed$new_chromosome)
#hg19_bed$new_chromosome[which(hg19_bed$new_chromosome=="Un")] <- 0


# read in + edit eqtl data

bg <- read_tsv("eqtl_data_metabrain/metabrain_basalganglia_eur_exposure_dat_snps_5kb_window.txt", col_types=cols(.default="c"))
bg$key_snp <- paste(bg$SNP_Chromosome,bg$SNP_Chromosome_Position,sep=":")
#bg$key_gene_start <- paste(bg$Gene_Chromosome,bg$start_position,sep=":")
#bg$key_gene_end <- paste(bg$Gene_Chromosome,bg$end_position,sep=":")

bg <- left_join(bg, setNames(hg19_bed[,c("key","new_chromosome","new_position")],c("key","hg19_snp_chr","hg19_snp_pos")), by=c("key_snp"="key"))
#bg <- left_join(bg, setNames(hg19_bed[,c("key","new_chromosome","new_position")],c("key","hg19_gene_chr","hg19_gene_start")), by=c("key_gene_start"="key"))
#bg <- left_join(bg, setNames(hg19_bed[,c("key","new_position")],c("key","hg19_gene_end")), by=c("key_gene_end"="key"))



cortex <- read_tsv("eqtl_data_metabrain/metabrain_cortex_eur_exposure_dat_snps_5kb_window.txt", col_types=cols(.default="c"))
cortex$key_snp <- paste(cortex$SNP_Chromosome,cortex$SNP_Chromosome_Position,sep=":")
#cortex$key_gene_start <- paste(cortex$Gene_Chromosome,cortex$start_position,sep=":")
#cortex$key_gene_end <- paste(cortex$Gene_Chromosome,cortex$end_position,sep=":")

cortex <- left_join(cortex, setNames(hg19_bed[,c("key","new_chromosome","new_position")],c("key","hg19_snp_chr","hg19_snp_pos")), by=c("key_snp"="key"))
#cortex <- left_join(cortex, setNames(hg19_bed[,c("key","new_chromosome","new_position")],c("key","hg19_gene_chr","hg19_gene_start")), by=c("key_gene_start"="key"))
#cortex <- left_join(cortex, setNames(hg19_bed[,c("key","new_position")],c("key","hg19_gene_end")), by=c("key_gene_end"="key"))

# write out

write.table(bg, "eqtl_data_metabrain/metabrain_basalganglia_eur_exposure_dat_snps_5kb_window_hg19.txt", sep = "\t", row.names = F, quote = F)
write.table(cortex, "eqtl_data_metabrain/metabrain_cortex_eur_exposure_dat_snps_5kb_window_hg19.txt", sep = "\t", row.names = F, quote = F)
