
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


### load exposure data

exp0 <- read_exposure_data(
  filename = "eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  phenotype_col = "GeneSymbol",
  samplesize_col = "NrSamples",
  min_pval = 1e-400
)

# keep only a subset of genes


    TO_REPLICATE <- exp0
    END <- length(unique(exp0$exposure))




exp_to_keep <- unique(exp0$exposure)
exp <- subset(exp0, (exp0$exposure %in% exp_to_keep) & !(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% TO_REPLICATE$exposure))


EXPOSURE_DAT <- gsub("_qc","", EXPOSURE_DATA)

eqtlgen <- read_tsv(gzfile("eqtl_data_eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"), col_types=cols(.default="c"))
eqtlgen$new_gene_id <- eqtlgen$Gene

eqtlgen_snp_counts <- data.frame(table(eqtlgen$SNP))
eqtlgen_snps_no_pleio <- subset(eqtlgen_snp_counts, eqtlgen_snp_counts$Freq == 1)

exp <- subset(exp, !(exp$SNP %in% eqtlgen_snps_no_pleio$Var1))


### load outcome data

#dat_keep <- readr::read_csv(str_c("full_results/significant_results_",OUTCOME,".txt"), col_types=cols(.default="c"))
dat_keep <- readr::read_csv(str_c("full_results/significant_results_replication_risk.txt"), col_types=cols(.default="c"))

dat_keep <- subset(dat_keep, dat_keep$tissue == EXPOSURE_DAT)

if (plyr::empty(dat_keep) == TRUE) {
    print("no significant genes in the discovery")
    print("mission_complete")
}


exp <- subset(exp, exp$exposure %in% dat_keep$exposure)


if (plyr::empty(exp) == TRUE) {
    print("exp is empty")
    print("mission_complete")
    exp <- data.frame()
}
