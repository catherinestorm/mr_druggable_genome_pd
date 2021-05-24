# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")
library("readr")

# environmental variables

EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME <- Sys.getenv("OUTCOME")

START <- 1
END <- NA


# problematic outcomes
REMOVE <- read.table("exposures_to_remove.txt", sep = ",", header = T, colClasses = "character")
REMOVE1 <- subset(REMOVE, REMOVE$outcome == OUTCOME & REMOVE$exposure_data == EXPOSURE_DATA)


# load exposure data

exp0 <- read_exposure_data(
  filename = "eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt",
  sep = "\t",
  snp_col = "Rsid",
  beta_col = "regression_slope",
  se_col = "se",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "nominal_pval",
  phenotype_col = "new_gene_id",
  min_pval = 1e-400
)

# add sample size

exp0$samplesize.exposure <- 1387


# keep only a subset of genes


TO_REPLICATE <- exp0
END <- length(unique(exp0$exposure))




exp_to_keep <- unique(exp0$exposure)
exp <- subset(exp0, (exp0$exposure %in% exp_to_keep) & !(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% TO_REPLICATE$exposure))



EXPOSURE_DAT <- gsub("_qc","", EXPOSURE_DATA)

psychencode <- read_tsv("eqtl_data_psychencode/DER-08a_hg19_eQTL.significant.txt", col_types=cols(.default="c"))
psychencode$new_gene_id <- gsub("\\..*", "", psychencode$gene_id)
alleles <- read_tsv("eqtl_data_psychencode/SNP_Information_Table_with_Alleles.txt", col_types=cols(.default="c"))
psychencode_full <- left_join(psychencode, alleles, by = c("SNP_id" = "PEC_id"))


psychencode_snp_counts <- data.frame(table(psychencode_full$Rsid))
psychencode_snps_no_pleio <- subset(psychencode_snp_counts, psychencode_snp_counts$Freq == 1)

exp <- subset(exp, !(exp$SNP %in% psychencode_snps_no_pleio$Var1))



#dat_keep <- readr::read_csv(str_c("full_results/significant_results_",OUTCOME,".txt"), col_types=cols(.default="c"))
dat_keep <- readr::read_csv(str_c("full_results/significant_results_replication_risk.txt"), col_types=cols(.default="c"))

dat_keep <- subset(dat_keep, dat_keep$tissue == EXPOSURE_DAT)

druggable <- read_tsv("druggable_genome_new.txt", col_types=cols(.default="c"))
dat_keep <- subset(druggable, druggable$gene_display_label %in% dat_keep$exposure)


if (plyr::empty(dat_keep) == TRUE) {
    print("no significant genes in the discovery")
    print("mission_complete")
    exp <- data.frame()
}

exp <- subset(exp, exp$exposure %in% dat_keep$gene_stable_id)

if (plyr::empty(exp) == TRUE) {
    print("exp is empty")
    print("mission_complete")
    exp <- data.frame()
}
