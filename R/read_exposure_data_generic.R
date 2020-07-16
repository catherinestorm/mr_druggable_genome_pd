# generic Rscript to read exposure data

## read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")


## environmental variables

EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME <- Sys.getenv("OUTCOME")

START <- as.numeric(Sys.getenv("START"))
END <- as.numeric(Sys.getenv("END"))


## problematic outcomes
REMOVE <- read.table("exposures_to_remove.txt", sep = ",", header = T, colClasses = "character")
REMOVE1 <- subset(REMOVE, REMOVE$outcome == OUTCOME & REMOVE$exposure_data == EXPOSURE_DATA)


## load exposure data

exp0 <- read_exposure_data(
  filename = "type in path to the file you with to use",
  sep = "separator between columns e.g. '/t' or ',' ",
  snp_col = "type in the name of the column that includes RS identifiers",
  beta_col = "type in the name of the column that includes the beta regression coefficients",
  se_col = "type in the name of the column that includes the standard errors",
  effect_allele_col = "type in the name of the column that includes the effect/alternate allele",
  other_allele_col = "type in the name of the column that includes the other/reference allele",
  pval_col = "type in the name of the column that includes the pvalue",
  phenotype_col = "type in the name of the column that includes the name of the gene",
  samplesize_col = "type in the name of the column that includes the number of samples used"
)


## keep only a subset of genes

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
