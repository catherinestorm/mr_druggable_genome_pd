### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                                     filename = "outcome_data/pd_replication_aao.txt",
                                     sep = "\t",
                                     snp_col = str_c("rsid_", EXPOSURE_DATA),
                                     beta_col = "Effect",
                                     se_col = "StdErr",
                                     effect_allele_col = "Allele1",
                                     other_allele_col = "Allele2",
                                     eaf_col = "Freq1",
                                     pval_col = "P-value")

out$outcome <- "pd_age_at_onset"

out$samplesize.outcome <- 17996
