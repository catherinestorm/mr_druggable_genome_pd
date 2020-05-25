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
  filename = "eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt",
  sep = "\t",
  snp_col = "Rsid",
  beta_col = "regression_slope",
  se_col = "se",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "nominal_pval",
  phenotype_col = "gene.exposure"
)

# add sample size

exp0$samplesize.exposure <- 1387


# keep only a subset of genes 

if (startsWith(OUTCOME, "replication_") == TRUE) {
    DISCOVERY_OUTCOME <- Sys.getenv("DISCOVERY_OUTCOME")
    TO_REPLICATE <- read.table(str_c(EXPOSURE_DATA, "_", DISCOVERY_OUTCOME, "/results/full_results_liberal_r2_0.2_", EXPOSURE_DATA,"_", DISCOVERY_OUTCOME, "_significant.txt"), sep = "\t", header = T, colClasses = "character")
    END <- length(unique(exp0$exposure))
} else if (is.na(END) == TRUE) {
    TO_REPLICATE <- exp0
    END <- length(unique(exp0$exposure))
} else {
    TO_REPLICATE <- exp0
}


exp_to_keep <- unique(exp0$exposure)[START:END]
exp <- subset(exp0, (exp0$exposure %in% exp_to_keep) & !(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% TO_REPLICATE$exposure))

