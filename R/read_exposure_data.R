# make the read_exposure data scripts

# if you have your own exposure data:
export FILENAME="FILENAME please include the file path"
export RSID="type in the name of the column that includes RS identifiers"
export BETA="type in the name of the column that inclufes the beta regression coefficients"
export STANDARD_ERROR="STANDARD_ERROR"
export EFFECT_ALLELE="EFFECT_ALLELE"
export OTHER_ALLELE="OTHER_ALLELE"
export PVALUE="PVALUE"
export GENE_NAME="GENE_NAME"
export SAMPLE_SIZE_COL="SAMPLE_SIZE_COL"
export SAMPLE_SIZE_DEFAULT="SAMPLE_SIZE_DEFAULT"
export SEPARATOR="SEPARATOR"
export EAF="EAF"



# read in packages

library("TwoSampleMR", lib.loc="/data/kronos/cstorm/mr")
library("dplyr")
library("stringr")


# environmental variables
FILENAME <- Sys.getenv("FILENAME")
RSID <- Sys.getenv("RSID")
BETA <- Sys.getenv("BETA")
STANDARD_ERROR <- Sys.getenv("STANDARD_ERROR")
EFFECT_ALLELE <- Sys.getenv("EFFECT_ALLELE")
OTHER_ALLELE <- Sys.getenv("OTHER_ALLELE")
PVALUE <- Sys.getenv("PVALUE")
GENE_NAME <- Sys.getenv("GENE_NAME")
SAMPLE_SIZE_COL <- Sys.getenv("SAMPLE_SIZE_COL")
SAMPLE_SIZE_DEFAULT <- as.numeric(Sys.getenv("SAMPLE_SIZE_DEFAULT"))
SEPARATOR <- Sys.getenv("SEPARATOR")
EAF <- Sys.getenv("EAF")


EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME <- Sys.getenv("OUTCOME")

START <- as.numeric(Sys.getenv("START"))
END <- as.numeric(Sys.getenv("END"))

# problematic outcomes
REMOVE <- read.table("exposures_to_remove.txt", sep = ",", header = T, colClasses = "character")
REMOVE1 <- subset(REMOVE, REMOVE$outcome == OUTCOME & REMOVE$exposure_data == EXPOSURE_DATA)



# load exposure data

exp0 <- read_exposure_data(
  filename = FILENAME,
  sep = SEPARATOR,
  snp_col = RSID,
  beta_col = BETA,
  se_col = STANDARD_ERROR,
  effect_allele_col = EFFECT_ALLELE,
  eaf_col = EAF,
  other_allele_col = OTHER_ALLELE,
  pval_col = PVALUE,
  phenotype_col = GENE_SAME,
  sample_size_col = SAMPLE_SIZE_COL
)

# add sample size

exp0$samplesize.exposure <- SAMPLE_SIZE_DEFAULT


# keep only a subset of genes 

exp_to_keep <- unique(exp0$exposure)[START:END]
exp <- subset(exp0, exp0$exposure %in% exp_to_keep)



