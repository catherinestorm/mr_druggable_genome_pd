# Put the Steiger-filtered MR input data into one file.


library(readr)
suppressMessages(library(dplyr, warn.conflict = FALSE, quietly = TRUE))

# find all Steiger-filtered MR input data generated
dat_all <- list.files(pattern = "steiger_stats", recursive = T, full.names = T)

# populate a data frame
dat_steiger_all <- data.frame()

for (i in 1:length(dat_all)) {
temp <- read_tsv(dat_all[i], col_types = cols())
temp$tissue <- gsub(".*liberal_r2_0.2_", "", dat_all[i])
temp$tissue <- gsub("\\_.*", "", temp$tissue)
temp$outcome_tissue <- paste(temp$outcome, temp$tissue, sep = "_")
temp$tissue_outcome <- paste(temp$tissue, temp$outcome, sep = "_")
dat_steiger_all <- distinct(rbind(dat_steiger_all, temp))
}

# populate a data frame
steiger_stats_summary <- as.data.frame.matrix(table(dat_steiger_all$outcome_tissue,dat_steiger_all$correct_causal_direction))

 steiger_stats_summary <- steiger_stats_summary[order(steiger_stats_summary[,1]),]
summary(steiger_stats_summary[,1])

steiger_stats_summary <- cbind(steiger_stats_summary, rowSums(steiger_stats_summary))
steiger_stats_summary <- cbind(row.names(steiger_stats_summary), steiger_stats_summary)
names(steiger_stats_summary) <- c("outcome","steiger_failed","steiger_passed","total_tests")
steiger_stats_summary$tissue <- gsub(".*\\_", "", steiger_stats_summary$outcome)
steiger_stats_summary$outcome <- gsub("_eqtlgen", "", steiger_stats_summary$outcome)

steiger_stats_summary <- steiger_stats_summary[c("outcome", "tissue", "steiger_failed", "steiger_passed", "total_tests")]

steiger_stats_summary$percent_failed <- round(steiger_stats_summary$steiger_failed / steiger_stats_summary$total_tests, digits=2)

# write it out
write.table(dat_steiger_all, "full_results/steiger_stats_liberal_all.txt", row.names = F, sep = "\t")
write.table(steiger_stats_summary, "full_results/steiger_stats_liberal_all_summary.txt", row.names = F, sep = "\t")
