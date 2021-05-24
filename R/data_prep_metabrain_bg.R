

library("dplyr")
library("readr")
library("biomaRt")

# load druggable genome

druggable <- read_tsv("druggable_genome_new.txt", col_types=cols(.default="c"))


# load basal ganglia data

bg <- read_csv("eqtl_data_metabrain/metabrain_basalganglia_eur.csv", col_types=cols(.default="c"))
bg$new_gene_id <- gsub("\\..*","", bg$Gene)

bg <- subset(bg, bg$new_gene_id %in% druggable$gene_stable_id)


bg$rsid <- gsub(".*:rs", "rs",bg$`SNP Name`)
bg$rsid <- gsub("\\:.*", "",bg$rsid)

bg$beta <- gsub("\\ .*", "",bg$`Meta-analysis Beta (s.e.)`)
bg$se <- gsub(".*\\(", "",bg$`Meta-analysis Beta (s.e.)`)
bg$se <- gsub(")", "",bg$se)

names(bg) <- gsub(" ", "_", names(bg))


# find gene start and end positions #grch38
biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)

my_ids <- data.frame(ensembl_gene_id_version=bg$Gene)
my_ids$ensembl_gene_id <- gsub("\\..*","", my_ids$ensembl_gene_id_version)

my_ids.version <- merge(my_ids, t2g, by= 'ensembl_gene_id')

bg_with_gene_pos <- left_join(bg, my_ids.version[,3:6], by = c("Gene"="ensembl_gene_id_version.y"))

bg_with_gene_pos <- bg_with_gene_pos[!is.na(bg_with_gene_pos$start_position),]



# keep 5kb window around the gene

if (all(bg_with_gene_pos$SNP_Chromosome == bg_with_gene_pos$Gene_Chromosome)) {
    bg_with_gene_pos_cis <- bg_with_gene_pos
} else {
    bg_with_gene_pos_cis <- bg_with_gene_pos[which(bg_with_gene_pos$SNP_Chromosome == bg_with_gene_pos$Gene_Chromosome),]
}

bg_with_gene_pos_cis$SNP_Chromosome_Position <- as.numeric(bg_with_gene_pos_cis$SNP_Chromosome_Position)

bg_5kb_window <- bg_with_gene_pos_cis[which(
    bg_with_gene_pos_cis$SNP_Chromosome_Position >= (bg_with_gene_pos_cis$start_position-5000) &
    bg_with_gene_pos_cis$SNP_Chromosome_Position <= (bg_with_gene_pos_cis$end_position+5000)
),]


bg_5kb_window <- tidyr::separate(bg_5kb_window,
    col="SNP_Alleles",
    "into"=c("allele1","allele2"),
    sep="/"
)
bg_5kb_window$other_allele <- NA
bg_5kb_window$other_allele[which(bg_5kb_window$allele1==bg_5kb_window$Effect_Allele)] <- bg_5kb_window$allele2[which(bg_5kb_window$allele1==bg_5kb_window$Effect_Allele)]
bg_5kb_window$other_allele[which(bg_5kb_window$allele2==bg_5kb_window$Effect_Allele)] <- bg_5kb_window$allele1[which(bg_5kb_window$allele2==bg_5kb_window$Effect_Allele)]

write.table(bg_5kb_window, "eqtl_data_metabrain/metabrain_basalganglia_eur_exposure_dat_snps_5kb_window.txt", sep = "\t", row.names = F, quote = F)


print("mission complete")
