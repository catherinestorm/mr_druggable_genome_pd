

# harmonise

dat0 <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out)

dat <- subset(dat0, dat0$mr_keep == TRUE)



# liberal analysis

dat0.2 <- clump_data(dat, clump_r2 = 0.2)


#STEIGER TEST


if (grepl("pd_risk_discovery", unique(dat0.2$outcome), fixed = TRUE) == TRUE) {
  dat0.2$r.outcome <- get_r_from_lor(dat0.2$beta.outcome, af = dat0.2$eaf.outcome, ncase = 13708, ncontrol = 95282, prevalence = 0.01)
} else if (grepl("pd_replication", unique(dat0.2$outcome), fixed = TRUE) == TRUE) {
  dat0.2$r.outcome <- get_r_from_lor(dat0.2$beta.outcome, af = dat0.2$eaf.outcome, ncase = 8036, ncontrol = 5803, prevalence = 0.01)
}



steiger_nalls0.2 <- directionality_test(dat0.2)

steiger_nalls0.2_1 <- subset(steiger_nalls0.2, correct_causal_direction == TRUE)

dat_steiger0.2 <- subset(dat0.2, dat0.2$exposure %in% steiger_nalls0.2_1$exposure)

write.table(steiger_nalls0.2, str_c(EXPOSURE_DATA, "_",OUTCOME, "/data/", "steiger_stats_liberal_r2_0.2_",EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), col.names=T, row.names = F, sep = "\t")

write.table(dat_steiger0.2, str_c(EXPOSURE_DATA, "_",OUTCOME, "/data/", "dat_steiger_liberal_r2_0.2_",EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), col.names=T, row.names = F, sep = "\t")



# MR analysis uncorrelated
two_snps_or_less <- as.data.frame(table(dat_steiger0.2$exposure))
dat_steiger0.2_two_snps_or_less <- subset(dat_steiger0.2, dat_steiger0.2$exposure %in% two_snps_or_less$Var1[which(two_snps_or_less$Freq <=2)])
mr_mrbase0.2 <- mr(dat_steiger0.2_two_snps_or_less)

if (plyr::empty(mr_mrbase0.2) == FALSE) {
  mr_mrbase0.2_keep <- mr_mrbase0.2[,c("exposure", "outcome", "nsnp", "method", "b", "se", "pval")]
  names(mr_mrbase0.2_keep) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
  mr_mrbase0.2_keep$clump_thresh <- "0.2"

  all_res <- mr_mrbase0.2_keep

  data_for_qval <- which(all_res$method == "IVW" | all_res$method == "Inverse variance weighted" | all_res$method == "Wald ratio")
  all_res[data_for_qval, "fdr_qval"] <- p.adjust(all_res[data_for_qval, "p"], method = "fdr")

  write.table(all_res, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), sep = "\t", row.names = F)
}




# correlated ivw where > 2 SNPs available

dat_steiger_keep <- subset(dat_steiger0.2, dat_steiger0.2$exposure %in% two_snps_or_less$Var1[which(two_snps_or_less$Freq > 2)])



if (plyr::empty(dat_steiger_keep) == TRUE) {

    print("no genes in this script have > 2 snps available")

} else {

unique_exposures0 <- unique(dat_steiger_keep$exposure)

unique_exposures0 <- sort(unique_exposures0)

unique_exposures <- unique_exposures0[!is.na(unique_exposures0)]


corr_results <- data.frame()

corr_egger_intercept <- data.frame()



for (i in 1:length(unique_exposures)) {
  dat_steiger_keep1 <- dat_steiger_keep[which(dat_steiger_keep$exposure == unique_exposures[i]), ] #subset data to keep only 1 exposure of interest

  dat2 <- dat_to_MRInput(dat_steiger_keep1, get_correlation=TRUE)

  ivw <- MendelianRandomization::mr_ivw(dat2[[1]], correl=TRUE)
  egger <- MendelianRandomization::mr_egger(dat2[[1]], correl=TRUE)
  maxlik <- MendelianRandomization::mr_maxlik(dat2[[1]], correl=TRUE)
  try(temp_data <- rbind(c(ivw@Exposure, ivw@Outcome, ivw@SNPs, ivw@class[1], ivw@Estimate, ivw@StdError, ivw@Pvalue),
                     c(egger@Exposure, egger@Outcome, egger@SNPs, egger@class[1], egger@Estimate, egger@StdError.Est, egger@Pvalue.Est),
                     c(maxlik@Exposure, maxlik@Outcome, maxlik@SNPs, maxlik@class[1], maxlik@Estimate, maxlik@StdError, maxlik@Pvalue)))

  try(corr_results <- rbind(corr_results, temp_data))

  try(temp_data1 <- cbind(egger@Exposure, egger@Outcome, egger@SNPs, egger@Intercept, egger@CILower.Int, egger@CIUpper.Int, egger@Pvalue.Int, egger@Heter.Stat[1], egger@Heter.Stat[2]))
  try(corr_egger_intercept <- rbind(corr_egger_intercept, temp_data1))


}



#correlated results
corr_results1 <- corr_results
names(corr_results1) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")

corr_results1[,5:7] <- sapply(corr_results1[,5:7], function(x){as.numeric(as.character(x))})


# combine results
corr_results2 <- corr_results1[,c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")]
corr_results2$clump_thresh <- "0.2"
all_res <- distinct(rbind(mr_mrbase0.2_keep, corr_results2))

write.table(all_res, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), sep = "\t", row.names = F)



#egger intercept and q values
corr_egger_intercept1 <- corr_egger_intercept
names(corr_egger_intercept1) <- c("exposure", "outcome", "nsnp", "intercept", "lower_ci", "upper_ci", "pvalue", "q", "q_pval")
corr_egger_intercept1[,3:9] <- sapply(corr_egger_intercept1[,3:9], function(x){as.numeric(as.character(x))})

cochrans_q <- corr_egger_intercept1

cochrans_q$nsnp <- as.numeric(as.character(cochrans_q$nsnp))
cochrans_q$q_df <- cochrans_q$nsnp -1

cochrans_q[which(cochrans_q$q >= cochrans_q$q_df), "i2"] <- (cochrans_q[which(cochrans_q$q >= cochrans_q$q_df), "q"] - cochrans_q[which(cochrans_q$q >= cochrans_q$q_df), "q_df"]) / cochrans_q[which(cochrans_q$q >= cochrans_q$q_df), "q"]
cochrans_q[which(cochrans_q$q < cochrans_q$q_df), "i2"] <- round(0, digits = 1)

cochrans_q$conf_int <- paste(round(cochrans_q$lower_ci, digits = 4), round(cochrans_q$upper_ci, digits = 4), sep = ", ")

qc_write_out <- cochrans_q[,c("exposure", "outcome", "nsnp", "intercept", "conf_int", "pvalue", "q", "q_pval", "i2")]
names(qc_write_out) <- c("exposure", "outcome", "nsnp", "egger_intercept", "egger_intercept_95_ci", "egger_intercept_pvalue", "cochrans_q", "cochrans_q_pval", "i2")

write.table(qc_write_out, str_c(EXPOSURE_DATA, "_",OUTCOME, "/results/", "results_liberal_r2_0.2_egger_cochransq_",EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), sep = "\t", row.names = F)

}

print("mission_complete")
