# generic Rscript to read outcome data

### read in outcome data

out <- read_outcome_data(snps = exp$SNP,
                                       filename = "type in path to the file you with to use,
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
