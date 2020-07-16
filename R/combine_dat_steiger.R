library(readr)
suppressMessages(library(dplyr, warn.conflict = FALSE, quietly = TRUE))

dat_all0 <- list.files(pattern = "dat_steiger_liberal", recursive = T, full.names = T)
dat_all <- dat_all0[780:790]

dat_steiger_all <- data.frame()

for (i in 1:length(dat_all)) {
temp <- read_tsv(dat_all[i], col_types = cols())
temp <- temp[,!(names(temp) == "chr.exposure")]
temp$tissue <- gsub(".*liberal_r2_0.2_", "", dat_all[i])
temp$tissue <- gsub("\\_.*", "", temp$tissue)
dat_steiger_all <- distinct(rbind(dat_steiger_all, temp))
}

write.table(dat_steiger_all, "full_results/dat_steiger_liberal_all.txt", row.names = F, sep = "\t")
