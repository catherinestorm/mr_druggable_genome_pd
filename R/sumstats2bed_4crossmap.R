#sumstats2bed_4crossmap.R

library("dplyr")
library("readr")

bg <- read_tsv("eqtl_data_metabrain/metabrain_basalganglia_eur_exposure_dat_snps_5kb_window.txt", col_types=cols(.default="c"))


cortex <- read_tsv("eqtl_data_metabrain/metabrain_cortex_eur_exposure_dat_snps_5kb_window.txt", col_types=cols(.default="c"))

bed <- setNames(cortex[,c("SNP_Chromosome","SNP_Chromosome_Position")], c("chr","pos")) %>%
#rbind(setNames(cortex[,c("Gene_Chromosome","start_position")], c("chr","pos"))) %>%
#rbind(setNames(cortex[,c("Gene_Chromosome","end_position")], c("chr","pos"))) %>%
rbind(setNames(bg[,c("SNP_Chromosome","SNP_Chromosome_Position")], c("chr","pos")))# %>%
#rbind(setNames(bg[,c("Gene_Chromosome","start_position")], c("chr","pos"))) %>%
#rbind(setNames(bg[,c("Gene_Chromosome","end_position")], c("chr","pos")))

bed <- distinct(bed)
bed$chr <- paste("chr",bed$chr,sep="")
bed$pos_minus1 <- as.integer(as.numeric(bed$pos)-1)
bed$pos_actual <- bed$pos

write.table(bed[,c("chr","pos_minus1","pos_actual")], "eqtl_data_metabrain/metabrain_4liftover_bed.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
