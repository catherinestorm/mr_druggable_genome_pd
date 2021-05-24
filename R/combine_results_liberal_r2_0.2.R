

library("stringr")
library("readr")
library("dplyr")

EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME <- Sys.getenv("OUTCOME")


###### read in results ######

file_list <- list.files(pattern = "results_liberal_r2_0.2_", full.names = T)
file_list_res <- str_subset(file_list, "^(?!.*egger)")
file_list_res <- str_subset(file_list_res, "^(?!.*full)")

# read in  each file and populate a data frame
full_results <- data.frame()

for (i in 1:length(file_list_res)) {
  temp_data <- read_tsv(file_list_res[i])
  temp_data <- temp_data[,!(names(temp_data)=="fdr_qval")]
  full_results <- rbind(full_results, temp_data)
}

# add FDR correction for the number of genes tested (applied to Wald ratio/IVW methods)
full_results <- as.data.frame(full_results)
full_results$fdr_qval <- NA
data_for_qval <- which(full_results$method == "IVW" | full_results$method == "Inverse variance weighted" | full_results$method == "Wald ratio")
full_results[data_for_qval, "fdr_qval"] <- p.adjust(full_results[data_for_qval,"p"], method = "fdr")


# for replication data, keep all genes that are nominally significant using the wald ratio/IVW
# for discovery data, keep all genes that are significant after multiple testing using the wald ratio/IVW
if (startsWith(OUTCOME, "replication_") == TRUE) {
wald_ivw_sign <- full_results[data_for_qval,]
wald_ivw_sign <- wald_ivw_sign[as.numeric(wald_ivw_sign$p) < 0.05,]
full_results_significant <- subset(full_results, full_results$exposure %in% wald_ivw_sign$exposure)
} else {
full_results_significant <- subset(full_results, as.numeric(full_results$fdr_qval) < 0.05)
}

# keep the results from all MR methods for the significant outcomes
full_results_significant <- full_results[full_results$exposure %in% full_results_significant$exposure, ]

# write out the data
write.table(full_results, str_c("full_results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, ".txt"), sep = "\t", row.names = F)

write.table(full_results_significant, str_c("full_results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_significant.txt"), sep = "\t", row.names = F)



###### read in QC ######

file_list_qc <- str_subset(file_list, "egger")
file_list_qc <- str_subset(file_list_qc, "^(?!.*full)")


if (length(file_list_qc) > 0) {
    full_results_qc <- data.frame()

    for (i in 1:length(file_list_qc)) {
      temp_data <- read_tsv(file_list_qc[i])
      full_results_qc <- rbind(full_results_qc, temp_data)
    }

    write.table(full_results_qc, str_c("full_results_liberal_r2_0.2_egger_cochransq_", EXPOSURE_DATA, "_", OUTCOME, ".txt"), sep = "\t", row.names = F)
}
# read in  each file and populate a data frame





paste(length(unique(full_results$exposure)), "genes tested")
paste(length(unique(full_results_significant$exposure)), "genes significant")
print("mission_complete")
