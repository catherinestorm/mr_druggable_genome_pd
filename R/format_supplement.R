

### this script will put the results and quality control metrics for all outcomes into one data frame

library(dplyr)
library(readr)


# FOR SIGNIFICANT GENES

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


# read in druggable genome data and add druggability tier to the data
druggable_new <- read.csv("druggable_genome_new.txt", sep = "\t", header = T, colClasses = "character")

significant_res_with_druggability_info <- left_join(significant_res,
                                         druggable_new[,c("gene_display_label","priority","chr_name", "gene_start", "gene_end")],
                                         by = c("exposure" = "gene_display_label"))


significant_res_with_druggability_info$tissue <- gsub("psychencode", "brain_eqtl", significant_res_with_druggability_info$tissue)
significant_res_with_druggability_info$tissue <- gsub("eqtlgen", "blood_eqtl", significant_res_with_druggability_info$tissue)
significant_res_with_druggability_info$tissue <- gsub("replication_pqtl", "blood_pqtl", significant_res_with_druggability_info$tissue)

names(significant_res_with_druggability_info)[names(significant_res_with_druggability_info) == "tissue"] <- "tissue_qtl_type"

write.table(significant_res_with_druggability_info, "full_results/significant_genes_results_all_outcomes.txt", row.names = F, sep = "\t")



# populate a data frame with the results for signficant genes across all outcomes
significant_qc <- data.frame()

for (i in 1:length(significant_qc_files)) {
  temp <- read_csv(significant_qc_files[i], col_types = cols())
  significant_qc <- distinct(rbind(significant_qc,temp))
}

significant_qc$tissue <- gsub("psychencode", "brain_eqtl", significant_qc$tissue)
significant_qc$tissue <- gsub("eqtlgen", "blood_eqtl", significant_qc$tissue)
significant_qc$tissue <- gsub("replication_pqtl", "blood_pqtl", significant_qc$tissue)

names(significant_qc)[names(significant_qc) == "tissue"] <- "tissue_qtl_type"

write.table(significant_qc, "full_results/significant_genes_qc_all_outcomes.txt", row.names = F, sep = "\t")









# FOR ALL RESULTS
# make a list of the files we want to combine
full_files <- list.files(path = "full_results", pattern = "full_results*",full.names = TRUE)
full_res_files <- full_files[!(grepl("qc", full_files))]
full_qc_files <- full_files[(grepl("qc", full_files))]



# populate a data frame with the results for signficant genes across all outcomes
full_res <- data.frame()

for (i in 1:length(full_res_files)) {
  temp <- read_csv(full_res_files[i], col_types = cols())
  full_res <- distinct(rbind(full_res,temp))
}

full_res <- subset(full_res, !(full_res$p == 0))

full_res$tissue <- gsub("psychencode", "brain_eqtl", full_res$tissue)
full_res$tissue <- gsub("eqtlgen", "blood_eqtl", full_res$tissue)
full_res$tissue <- gsub("replication_pqtl", "blood_pqtl", full_res$tissue)

names(full_res)[names(full_res) == "tissue"] <- "tissue_qtl_type"

write.table(full_res, "full_results/all_genes_results_all_outcomes.txt", row.names = F, sep = "\t")



# populate a data frame with the results for signficant genes across all outcomes
full_qc <- data.frame()

for (i in 1:length(full_qc_files)) {
  temp <- read_csv(full_qc_files[i], col_types = cols())
  full_qc <- distinct(rbind(full_qc,temp))
}


full_qc$tissue <- gsub("psychencode", "brain_eqtl", full_qc$tissue)
full_qc$tissue <- gsub("eqtlgen", "blood_eqtl", full_qc$tissue)
full_qc$tissue <- gsub("replication_pqtl", "blood_pqtl", full_qc$tissue)

names(full_qc)[names(full_qc) == "tissue"] <- "tissue_qtl_type"


write.table(full_qc, "full_results/all_genes_qc_all_outcomes.txt", row.names = F, sep = "\t")








print("mission complete")
