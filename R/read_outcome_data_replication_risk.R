### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                                     filename = "outcome_data/nalls2019_replication_risk.txt",
                                     sep = "\t",
                                     snp_col = str_c("rsid_", EXPOSURE_DATA),
                                     beta_col = "beta",
                                     se_col = "se",
                                     effect_allele_col = "Allele1",
                                     other_allele_col = "Allele2",
                                     eaf_col = "Freq1",
                                     pval_col = "P-value",
                                     samplesize_col = "Weight")
  
out$outcome <- "nalls2019"

#out$samplesize.outcome <- 13839