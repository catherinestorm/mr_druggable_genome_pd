library("readr")
library("plyr")
suppressMessages(library(dplyr, warn.conflict = FALSE, quietly = TRUE))
library("stringr")



EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME <- Sys.getenv("OUTCOME")

# read in files to edit
total <- read_csv(str_c("full_results/full_results_",OUTCOME,".txt"), col_types = cols())
sign_total  <- read_csv(str_c("full_results/significant_results_",OUTCOME,".txt"), col_types = cols())

total_qc <- read_csv(str_c("full_results/full_results_",OUTCOME,"_qc.txt"), col_types = cols())
sign_total_qc <- read_csv(str_c("full_results/significant_results_",OUTCOME,"_qc.txt"), col_types = cols())

final_results_report <- read_csv("full_results/final_results_report.txt", col_types = cols())


# read in druggable genome and keep the columns we need
druggable_genome <- read_tsv("druggable_genome_new.txt", col_types = cols())
druggable_genome <- druggable_genome[,c("gene_display_label", "priority")]
names(druggable_genome) <- c("exposure", "druggability_tier")





# read in results
exposure <- read_tsv(str_c(EXPOSURE_DATA, "_",OUTCOME ,"/results/full_results_",EXPOSURE_DATA,"_",OUTCOME,".txt"), col_types = cols())


# add druggability info & calculate odds ratio
exposure <- left_join(exposure, druggable_genome, by = "exposure")

exposure$or <- exp(exposure$beta)
exposure$or_lci95 <- exp(exposure$beta-1.96*exposure$se)
exposure$or_uci95 <- exp(exposure$beta+1.96*exposure$se)

exposure$beta_lci95 <- exposure$beta-1.96*exposure$se
exposure$beta_uci95 <- exposure$beta+1.96*exposure$se

total <- distinct(rbind(total, exposure))

write.table(total, str_c("full_results/full_results_",OUTCOME,".txt"), sep = ",", row.names = F)


### significant results
sign <- read_tsv(str_c(EXPOSURE_DATA, "_",OUTCOME ,"/results/full_results_",EXPOSURE_DATA,"_",OUTCOME,"_significant.txt"), col_types = cols())

# add druggability info & calculate odds ratio
if (empty(sign) == F) {

sign <- left_join(sign, druggable_genome, by = "exposure")

sign$or <- exp(sign$beta)
sign$or_lci95 <- exp(sign$beta-1.96*sign$se)
sign$or_uci95 <- exp(sign$beta+1.96*sign$se)

sign$beta_lci95 <- sign$beta-1.96*sign$se
sign$beta_uci95 <- sign$beta+1.96*sign$se

sign_total <- distinct(rbind(sign_total, sign))

write.table(sign_total, str_c("full_results/significant_results_",OUTCOME,".txt"), sep = ",", row.names = F)


}




# qc
exposure_qc <- read_tsv(str_c(EXPOSURE_DATA,"_",OUTCOME ,"/results/full_results_liberal_r2_0.2_egger_cochransq_",EXPOSURE_DATA,"_",OUTCOME,".txt"), col_types = cols())

total_qc <- distinct(rbind(total_qc, exposure_qc))
sign_total_qc_new <- subset(total_qc, total_qc$exposure %in% sign$exposure)
sign_total_qc <- distinct(rbind(sign_total_qc, sign_total_qc_new))

write.table(total_qc, str_c("full_results/full_results_",OUTCOME,"_qc.txt"), sep = ",", row.names = F)
write.table(sign_total_qc, str_c("full_results/significant_results_",OUTCOME,"_qc.txt"), sep = ",", row.names = F)




final_results_report_new <- data.frame(EXPOSURE_DATA, OUTCOME, length(unique(exposure$exposure)), length(unique(sign$exposure)))
names(final_results_report_new) <- names(final_results_report)

final_results_report <- distinct(rbind(final_results_report, final_results_report_new))
tot_count  <- data.frame("total", OUTCOME, length(unique(total$exposure)), length(unique(sign_total$exposure)))
names(tot_count) <- names(final_results_report)



if (empty(final_results_report[final_results_report$exposure == "total" & final_results_report$outcome == OUTCOME,]) == T) {
final_results_report <- rbind(final_results_report,tot_count)
} else {
final_results_report[final_results_report$exposure == "total" & final_results_report$outcome == OUTCOME,"n_tested"] <- length(unique(total$exposure))
final_results_report[final_results_report$exposure == "total" & final_results_report$outcome == OUTCOME,"n_significant"] <- length(unique(sign_total$exposure))
}

write.table(final_results_report, "full_results/final_results_report.txt", sep = ",", row.names = F)



