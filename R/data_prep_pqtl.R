# pqtl data prep for druggable genome replication


library(dplyr)
library(stringr)

# genes to replicate
druggable_genome <- read.csv("druggable_genome_new.txt", sep = "\t", header = T, stringsAsFactors = F)
sign <- as.data.frame(read_tsv("full_results/significant_genes_results_all_outcomes.txt"))
sign <- subset(sign, sign$outcome != "nalls2014")

#druggable_genome_sign <- druggable_genome
druggable_genome_sign <- subset(druggable_genome, druggable_genome$gene_display_label %in% sign$exposure)

length(unique(druggable_genome_sign$gene_display_label))


# read in pqtl data
# obtained from supplementary material of the relevant papers
# please see the original publications for these papers

# Emilsson et al. 2018. Science. Table S1, edited with alleles obtained from the authors.
emilsson2018 <- read.csv("pqtl_data/emilsson_2018_for_mr.csv", header = T, stringsAsFactors = F)
emilsson2018_replicate <- subset(emilsson2018, emilsson2018$Protein %in% druggable_genome_sign$gene_display_label)
emilsson2018_replicate$gene_and_source <- str_c(emilsson2018_replicate$Protein, "_emilsson2018")
emilsson2018_replicate_keep <- emilsson2018_replicate[,c("Lead.pSNP", "chr","pos","effect_allele", "other_allele","eaf","Beta.coeff","Stderr","P.value","Protein","gene_and_source","N.subjects")]
names(emilsson2018_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

# Suhre et al. 2017. Nat Comm. Supplementary Dataset 1. KORA
suhre_2017 <- read.csv("pqtl_data/suhre_2017.csv", header = T, stringsAsFactors = F)
suhre_2017_replicate <- subset(suhre_2017, suhre_2017$UniProt %in% druggable_genome_sign$accession)
suhre_2017_replicate$gene_and_source <- str_c(suhre_2017_replicate$EntrezGeneSymbol, "_suhre2017")
suhre_2017_replicate$position <- gsub(",","",suhre_2017_replicate$Position)
suhre_2017_replicate$maf <- gsub("%", "", suhre_2017_replicate$MAF)
suhre_2017_replicate$maf <- as.numeric(suhre_2017_replicate$maf)/100
suhre_2017_replicate_keep <- suhre_2017_replicate[,c("SNP", "Chr","position","Minor.Allele", "Major.Allele","maf","Beta.inv","S.E.","P.Value..inv.","EntrezGeneSymbol","gene_and_source","N")]
names(suhre_2017_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

# Sun et al. 2018. Nature. Supplementary Table 4.
sun_2018 <- read.csv("pqtl_data/sun_2018.csv", header = T, stringsAsFactors = F)
sun_2018_replicate <- subset(sun_2018, sun_2018$UniProt %in% druggable_genome_sign$accession)
sun_2018_replicate$Target[sun_2018_replicate$Target == "Cathepsin B"] <- "CTSB"
sun_2018_replicate$Target[sun_2018_replicate$Target == "DHPR"] <- "QDPR"
sun_2018_replicate$gene_and_source <- str_c(sun_2018_replicate$Target, "_sun2018")
sun_2018_replicate$position <- gsub(",","",sun_2018_replicate$Pos)
sun_2018_replicate$n <- 3301
sun_2018_replicate_keep <- sun_2018_replicate[,c("Sentinel.variant.", "Chr","position","Effect.Allele..EA.", "Other.Allele..OA.","EAF","meta_beta","meta_se","meta_p","Target","gene_and_source","n")]
names(sun_2018_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")




# combine
complete_pqtl_data_for_druggable_genome_replication <- emilsson2018_replicate_keep %>% rbind(hillary_2019_replicate_keep) %>% rbind(suhre_2017_replicate_keep) %>% rbind(sun_2018_replicate_keep)

complete_pqtl_data_for_druggable_genome_replication <- distinct(complete_pqtl_data_for_druggable_genome_replication)


# add gene info from druggable genome file
gene_chr_pos <- druggable_genome_sign[,c("gene_display_label","chr_name","gene_start","gene_end")]
names(gene_chr_pos) <- c("protein","gene_chr","gene_start","gene_end")

complete_pqtl_data_for_druggable_genome_replication1 <- left_join(complete_pqtl_data_for_druggable_genome_replication, gene_chr_pos, by = "protein")

# clarify if pQTL is acting in cis or trans
complete_pqtl_data_for_druggable_genome_replication1$cis_trans <- NA
complete_pqtl_data_for_druggable_genome_replication1$cis_trans[complete_pqtl_data_for_druggable_genome_replication1$chr == complete_pqtl_data_for_druggable_genome_replication1$gene_chr] <- "cis"
complete_pqtl_data_for_druggable_genome_replication1$cis_trans[complete_pqtl_data_for_druggable_genome_replication1$chr != complete_pqtl_data_for_druggable_genome_replication1$gene_chr] <- "trans"

# check how many proteins were covered
length(unique(complete_pqtl_data_for_druggable_genome_replication1$gene_and_source))

write.table(complete_pqtl_data_for_druggable_genome_replication, "pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt",sep ="\t",row.names = F)
