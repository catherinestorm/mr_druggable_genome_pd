
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

if (startsWith(OUTCOME, "replication_") == TRUE) {
    discovery_outcomes <- read.table("discovery_outcomes.txt")

    TO_REPLICATE <- data.frame()
    for (i in 1:nrow(discovery_outcomes)) {

        DISCOVERY_OUTCOME <- discovery_outcomes[i,]

        temp <- read.table(str_c(EXPOSURE_DATA, "_", DISCOVERY_OUTCOME, "/results/full_results_liberal_r2_0.2_", EXPOSURE_DATA,"_", DISCOVERY_OUTCOME, "_significant.txt"), sep = "\t", header = T, colClasses = "character")

        TO_REPLICATE <- distinct(rbind(TO_REPLICATE, temp))

    }

    END <- length(unique(exp0$exposure))
} else if (is.na(END) == TRUE) {
    TO_REPLICATE <- exp0
    END <- length(unique(exp0$exposure))
} else {
    TO_REPLICATE <- exp0
}



exp_to_keep <- unique(exp0$exposure)[START:END]
exp <- subset(exp0, (exp0$exposure %in% exp_to_keep) & !(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% TO_REPLICATE$exposure))
