
# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")


# environmental variables

EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME <- Sys.getenv("OUTCOME")

START <- as.numeric(Sys.getenv("START"))
END <- as.numeric(Sys.getenv("END"))


# problematic outcomes
REMOVE <- read.table("exposures_to_remove.txt", sep = ",", header = T, colClasses = "character")
REMOVE1 <- subset(REMOVE, REMOVE$outcome == OUTCOME & REMOVE$exposure_data == EXPOSURE_DATA)


# load exposure data

exp0 <- read_exposure_data(
  filename = "pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt",
  sep = "\t",
  snp_col = "snp",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "gene_and_source",
  samplesize_col = "sample_size",
  min_pval = 1e-400
)

## keep only a subset of genes

if (startsWith(OUTCOME, "replication_") == TRUE) {

    # for a replication outcome, we will read in the significant results from the discovery phase

    DISCOVERY_OUTCOME <- Sys.getenv("DISCOVERY_OUTCOME")
    TO_REPLICATE <- read.table(str_c(EXPOSURE_DATA, "_", DISCOVERY_OUTCOME, "/results/full_results_liberal_r2_0.2_", EXPOSURE_DATA,"_", DISCOVERY_OUTCOME, "_significant.txt"), sep = "\t", header = T, colClasses = "character")

    # subset the exposure data, keeping only the genes that reached significance in the discovery phase (clumping at 0.2)
    exp <- subset(exp0, (!(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% TO_REPLICATE$exposure)))

} else if ((startsWith(EXPOSURE_DATA, "replication_") == TRUE)) {

    exp_to_keep <- distinct(data.frame("exp_source" = exp0$exposure, "exp" = gsub("\\_.*","", exp0$exposure), stringsAsFactors = F))

    TO_REPLICATE <- read.table(str_c("full_results/significant_results_",OUTCOME,".txt"), sep = ",", header = T, colClasses = "character")

    # subset the exposure data, keeping only the genes that reached significance in the discovery phase (clumping at 0.2)

    exp_to_keep <- subset(exp_to_keep, exp_to_keep$exp %in% TO_REPLICATE$exposure)

    exp <- subset(exp0, (!(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% exp_to_keep$exp_source)))

    exp_to_keep <- exp_to_keep$exp_source

} else if (is.na(END) == TRUE) {

    # when clumping at 0.001, we will read in the significant results from the main mr analysis (clumping at 0.2). we will subset the exposure data later, depending on if this file is empty or not

    TO_REPLICATE <- read.table(str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/full_results_liberal_r2_0.2_", EXPOSURE_DATA,"_", OUTCOME, "_significant.txt"), sep = "\t", header = T, colClasses = "character")


} else {

    # when clumping at 0.2 for a discovery outcome which runs several scripts in parallel, we want to subset the exposure data to run only a chunk of the genes

    exp_to_keep <- unique(exp0$exposure)[START:END]
    exp <- subset(exp0, !(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% exp_to_keep))
}

if (plyr::empty(exp) == TRUE) {
    print("no pqtls available for genes in this outcome")
    print("mission_complete")
}
