

### this script will put the results and quality control metrics for all significant outcomes into one data frame

library(dplyr)
library(readr)



# make a list of the files we want to combine
significant_files <- list.files(path = "full_results", pattern = "significant_results*",full.names = TRUE)
significant_res_files <- significant_files[!(grepl("qc", significant_files))]
significant_qc_files <- significant_files[(grepl("qc", significant_files))]



# populate a data frame with the results for signficant genes across all outcomes
significant_res <- data.frame()

for (i in 1:length(significant_res_files)) {
  temp <- read_csv(significant_res_files[i], col_types = cols())
  significant_res <- distinct(rbind(significant_res,temp))
}

significant_res <- subset(significant_res, !(significant_res$p == 0))



druggable_new <- read.csv("druggable_genome_new.txt", sep = "\t", header = T, colClasses = "character")

significant_res_with_druggability_info <- left_join(significant_res, 
                                         druggable_new[,c("gene_display_label","priority","chr_name", "gene_start", "gene_end")],
                                         by = c("exposure" = "gene_display_label"))


write.table(significant_res_with_druggability_info, "full_results/significant_genes_results_all_outcomes.txt", row.names = F, sep = "\t")



# populate a data frame with the results for signficant genes across all outcomes
significant_qc <- data.frame()

for (i in 1:length(significant_qc_files)) {
  temp <- read_csv(significant_qc_files[i], col_types = cols())
  significant_qc <- distinct(rbind(significant_qc,temp))
}

write.table(significant_qc, "full_results/significant_genes_qc_all_outcomes.txt", row.names = F, sep = "\t")



print("mission complete")