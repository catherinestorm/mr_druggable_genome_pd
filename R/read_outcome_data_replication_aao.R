### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                                     filename = "outcome_data/blauwendraat2019_replication_aao.txt",
                                     sep = "\t",
                                     snp_col = str_c("rsid_", EXPOSURE_DATA),
                                     beta_col = "Effect",
                                     se_col = "StdErr",
                                     effect_allele_col = "Allele1",
                                     other_allele_col = "Allele2",
                                     eaf_col = "Freq1",
                                     pval_col = "P-value")
  
out$outcome <- "blauwendraat2019"

out$samplesize.outcome <- 17996
