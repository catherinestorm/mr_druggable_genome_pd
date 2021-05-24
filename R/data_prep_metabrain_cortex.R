

library("dplyr")
library("readr")
library("biomaRt")


# load druggable genome

druggable <- read_tsv("druggable_genome_new.txt", col_types=cols(.default="c"))


# load cortex data

cortex <- read_csv("eqtl_data_metabrain/metabrain_cortex_eur.csv", col_types=cols(.default="c"))
cortex$new_gene_id <- gsub("\\..*","", cortex$Gene)

cortex <- subset(cortex, cortex$new_gene_id %in% druggable$gene_stable_id)


cortex$rsid <- gsub(".*:rs", "rs",cortex$`SNP Name`)
cortex$rsid <- gsub("\\:.*", "",cortex$rsid)

cortex$beta <- gsub("\\ .*", "",cortex$`Meta-analysis Beta (s.e.)`)
cortex$se <- gsub(".*\\(", "",cortex$`Meta-analysis Beta (s.e.)`)
cortex$se <- gsub(")", "",cortex$se)

names(cortex) <- gsub(" ", "_", names(cortex))


# find gene start and end positions #grch38
biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)

my_ids <- data.frame(ensembl_gene_id_version=cortex$Gene)
my_ids$ensembl_gene_id <- gsub("\\..*","", my_ids$ensembl_gene_id_version)

my_ids.version <- merge(my_ids, t2g, by= 'ensembl_gene_id')

cortex_with_gene_pos <- left_join(cortex, my_ids.version[,3:6], by = c("Gene"="ensembl_gene_id_version.y"))

cortex_with_gene_pos <- cortex_with_gene_pos[!is.na(cortex_with_gene_pos$start_position),]


# keep 5kb window around the gene

if (all(cortex_with_gene_pos$SNP_Chromosome == cortex_with_gene_pos$Gene_Chromosome)) {
    cortex_with_gene_pos_cis <- cortex_with_gene_pos
} else {
    cortex_with_gene_pos_cis <- cortex_with_gene_pos[which(cortex_with_gene_pos$SNP_Chromosome == cortex_with_gene_pos$Gene_Chromosome),]
}

cortex_with_gene_pos_cis$SNP_Chromosome_Position <- as.numeric(cortex_with_gene_pos_cis$SNP_Chromosome_Position)

cortex_5kb_window <- cortex_with_gene_pos_cis[which(
    cortex_with_gene_pos_cis$SNP_Chromosome_Position >= (cortex_with_gene_pos_cis$start_position-5000) &
    cortex_with_gene_pos_cis$SNP_Chromosome_Position <= (cortex_with_gene_pos_cis$end_position+5000)
),]

cortex_5kb_window <- tidyr::separate(cortex_5kb_window,
    col="SNP_Alleles",
    "into"=c("allele1","allele2"),
    sep="/"
)
cortex_5kb_window$other_allele <- NA
cortex_5kb_window$other_allele[which(cortex_5kb_window$allele1==cortex_5kb_window$Effect_Allele)] <- cortex_5kb_window$allele2[which(cortex_5kb_window$allele1==cortex_5kb_window$Effect_Allele)]
cortex_5kb_window$other_allele[which(cortex_5kb_window$allele2==cortex_5kb_window$Effect_Allele)] <- cortex_5kb_window$allele1[which(cortex_5kb_window$allele2==cortex_5kb_window$Effect_Allele)]

write.table(cortex_5kb_window, "eqtl_data_metabrain/metabrain_cortex_eur_exposure_dat_snps_5kb_window.txt", sep = "\t", row.names = F, quote = F)


print("mission complete")
