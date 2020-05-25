

# harmonise

dat0 <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out)

dat <- subset(dat0, dat0$mr_keep == TRUE)



# liberal analysis

dat0.2 <- clump_data(dat, clump_r2 = 0.2)


#STEIGER TEST

steiger_nalls0.2 <- directionality_test(dat0.2)

steiger_nalls0.2_1 <- subset(steiger_nalls0.2, correct_causal_direction == TRUE)
  
dat_steiger0.2 <- subset(dat0.2, dat0.2$exposure %in% steiger_nalls0.2_1$exposure)

write.table(dat_steiger0.2, str_c(EXPOSURE_DATA, "_",OUTCOME, "/data/", "dat_steiger_liberal_r2_0.2_",EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), col.names=T, row.names = F, sep = "\t")

# MR analysis uncorrelated
mr_mrbase0.2 <- mr(dat_steiger0.2)
mr_mrbase0.2_keep <- subset(mr_mrbase0.2, mr_mrbase0.2$nsnp <= 2)

mr_mrbase0.2_keep <- mr_mrbase0.2_keep[,c("exposure", "outcome", "nsnp", "method", "b", "se", "pval")]
names(mr_mrbase0.2_keep) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")

write.table(mr_mrbase0.2_keep, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), sep = "\t", row.names = F)







# correlated where > 2 SNPs available

dat_steiger_keep <- subset(dat_steiger0.2, dat_steiger0.2$exposure %in% mr_mrbase0.2[which(mr_mrbase0.2$nsnp > 2), "exposure"])



if (plyr::empty(dat_steiger_keep) == TRUE) {

    print("no more genes in this script have > 2 snps available")

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
  temp_data <- rbind(c(ivw@Outcome, ivw@SNPs, ivw@class[1], ivw@Estimate, ivw@StdError, ivw@Pvalue),
                     c(egger@Outcome, egger@SNPs, egger@class[1], egger@Estimate, egger@StdError.Est, egger@Pvalue.Est),
                     c(maxlik@Outcome, maxlik@SNPs, maxlik@class[1], maxlik@Estimate, maxlik@StdError, maxlik@Pvalue))
  
  corr_results <- rbind(corr_results, temp_data)
  
  temp_data1 <- cbind(egger@Outcome, egger@Intercept, egger@CILower.Int, egger@CIUpper.Int, egger@Pvalue.Int, egger@Heter.Stat[1], egger@Heter.Stat[2])
  corr_egger_intercept <- rbind(corr_egger_intercept, temp_data1)
  
  
}



#correlated results
corr_results1 <- corr_results
names(corr_results1) <- c("outcome", "nsnp", "method", "beta", "se", "p")

corr_results1[which(corr_results1$method == "IVW"), "exposure"] <- unique_exposures
corr_results1[which(corr_results1$method == "Egger"), "exposure"] <- unique_exposures
corr_results1[which(corr_results1$method == "MaxLik"), "exposure"] <- unique_exposures
corr_results1[,4:6] <- sapply(corr_results1[,4:6], function(x){as.numeric(as.character(x))})


# combine results
corr_results2 <- corr_results1[,c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")]
all_res <- distinct(rbind(mr_mrbase0.2_keep, corr_results2))

write.table(all_res, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), sep = "\t", row.names = F)



#egger intercept and q values
corr_egger_intercept1 <- corr_egger_intercept
names(corr_egger_intercept1) <- c("Outcome", "Intercept", "lower_ci", "upper_ci", "pvalue", "Q", "Q_pval")
corr_egger_intercept1$exposure <- unique_exposures
corr_egger_intercept1[,2:7] <- sapply(corr_egger_intercept1[,2:7], function(x){as.numeric(as.character(x))})

cochrans_q <- merge(corr_egger_intercept1,
                    distinct(corr_results1[, c("exposure", "nsnp")]),
                    by = "exposure")

cochrans_q$nsnp <- as.numeric((cochrans_q$nsnp))
cochrans_q$Q_df <- cochrans_q$nsnp -1

cochrans_q[which(cochrans_q$Q >= cochrans_q$Q_df), "I2"] <- (cochrans_q[which(cochrans_q$Q >= cochrans_q$Q_df), "Q"] - cochrans_q[which(cochrans_q$Q >= cochrans_q$Q_df), "Q_df"]) / cochrans_q[which(cochrans_q$Q >= cochrans_q$Q_df), "Q"]
cochrans_q[which(cochrans_q$Q < cochrans_q$Q_df), "I2"] <- round(0, digits = 1)
cochrans_q$conf_int <- paste(round(cochrans_q$lower_ci, digits = 4), round(cochrans_q$upper_ci, digits = 4), sep = ", ")

qc_write_out <- cochrans_q[,c("exposure", "Outcome", "nsnp", "Intercept", "conf_int", "pvalue", "Q", "Q_pval", "I2")]
names(qc_write_out) <- c("exposure", "outcome", "nsnp", "egger_intercept", "egger_intercept_95_ci", "egger_intercept_pvalue", "cochrans_q", "cochrans_q_pval", "i2")

write.table(qc_write_out, str_c(EXPOSURE_DATA, "_",OUTCOME, "/results/", "results_liberal_r2_0.2_egger_cochransq_",EXPOSURE_DATA, "_", OUTCOME, "_", START, ".txt"), sep = "\t", row.names = F)

}

print("mission_complete")