# Put the Steiger-filtered MR input data into one file.


library(readr)
suppressMessages(library(dplyr, warn.conflict = FALSE, quietly = TRUE))

# find all Steiger-filtered MR input data generated
dat_all0 <- list.files(pattern = "dat_steiger_liberal", recursive = T, full.names = T)

# populate a data frame
dat_steiger_all <- data.frame()

for (i in 1:length(dat_all)) {
temp <- read_tsv(dat_all[i], col_types = cols())
temp <- temp[,!(names(temp) == "chr.exposure")]
temp$tissue <- gsub(".*liberal_r2_0.2_", "", dat_all[i])
temp$tissue <- gsub("\\_.*", "", temp$tissue)
temp <- temp[,c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure", "beta.outcome","eaf.exposure","eaf.outcome", "remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","pval.exposure","samplesize.exposure","se.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure","data_source.exposure","action","mr_keep","tissue")]
dat_steiger_all <- distinct(rbind(dat_steiger_all, temp))
}

# write it out
write.table(dat_steiger_all, "full_results/dat_steiger_liberal_all.txt", row.names = F, sep = "\t")
