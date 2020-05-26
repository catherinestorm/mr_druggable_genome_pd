# keep only the genes that were significant in the main mr analysis (clumping at 0.2)
TO_REPLICATE <- read.table(str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/full_results_liberal_r2_0.2_", EXPOSURE_DATA,"_", OUTCOME, "_significant.txt"), sep = "\t", header = T, colClasses = "character")


if (plyr::empty(TO_REPLICATE) == TRUE) {

liberal_full <- readr::read_tsv(str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, ".txt"))
liberal_sign <- readr::read_tsv(str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_significant.txt"))

liberal_full$clump_tresh <- "0.2"
liberal_sign$clump_tresh <- "0.2"

liberal_full$tissue <- EXPOSURE_DATA
liberal_sign$tissue <- EXPOSURE_DATA

write.table(liberal_full, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_", EXPOSURE_DATA, "_", OUTCOME, ".txt"), sep = "\t", row.names = F)

write.table(liberal_sign, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_", EXPOSURE_DATA, "_", OUTCOME, "_significant.txt"), sep = "\t", row.names = F)


message1 <- str_c(length(unique(liberal_full$exposure)), " genes tested in main analysis")
message2 <- str_c("no genes were significant using the main analysis (r2 = 0.2)")




} else {


exp <- subset(exp0, (exp0$exposure %in% exp_to_keep) & !(exp0$exposure %in% REMOVE1$exposure) & (exp0$exposure %in% TO_REPLICATE$exposure))



# harmonise

dat0 <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out)



# how many snps left?

dat <- subset(dat0, dat0$mr_keep == TRUE)

# strict analysis

dat0.001 <- clump_data(dat)


# STEIGER TEST

steiger_nalls0.001 <- directionality_test(dat0.001)

steiger_nalls0.001_1 <- subset(steiger_nalls0.001, correct_causal_direction == TRUE)
  
dat_steiger0.001 <- subset(dat0.001, dat0.001$exposure %in% steiger_nalls0.001_1$exposure)

write.table(dat_steiger0.001, str_c(EXPOSURE_DATA, "_",OUTCOME, "/data/", "dat_steiger_conservative_r2_0.001_",EXPOSURE_DATA, "_", OUTCOME, ".txt"), col.names=T, row.names = F, sep = "\t")



# MR analysis uncorrelated

mr_mrbase0.001 <- mr(dat_steiger0.001)
mr_mrbase0.001 <- mr_mrbase0.001[,c("exposure", "outcome", "nsnp", "method", "b", "se", "pval")]
names(mr_mrbase0.001) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")

# multiple testing correction
mr_mrbase0.001$fdr_qval <- "NA"
data_for_qval <- which(mr_mrbase0.001$method == "IVW" | mr_mrbase0.001$method == "Inverse variance weighted" | mr_mrbase0.001$method == "Wald ratio")
mr_mrbase0.001[data_for_qval, "fdr_qval"] <- p.adjust(mr_mrbase0.001[data_for_qval,"p"], method = "fdr")
mr_mrbase0.001$clump_tresh <- "0.001"



#if (startsWith(OUTCOME, "replication_") == TRUE) {
mr_mrbase0.001_sign <- mr_mrbase0.001[mr_mrbase0.001$p < 0.05,]
#} else {
#mr_mrbase0.001_sign <- mr_mrbase0.001[mr_mrbase0.001$fdr_qval < 0.05,]
#}

write.table(mr_mrbase0.001, str_c(EXPOSURE_DATA, "_",OUTCOME, "/results/", "full_results_conservative_r2_0.001_",EXPOSURE_DATA, "_", OUTCOME, ".txt"), col.names=T, row.names = F, sep = "\t")


# combine with other data
liberal_full <- readr::read_tsv(str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, ".txt"))

liberal_full$clump_tresh <- "0.2"
liberal_full <- distinct(rbind(mr_mrbase0.001,liberal_full))
liberal_full$tissue <- EXPOSURE_DATA
write.table(liberal_full, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_", EXPOSURE_DATA, "_", OUTCOME, ".txt"), sep = "\t", row.names = F)


liberal_sign0 <- readr::read_tsv(str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_liberal_r2_0.2_", EXPOSURE_DATA, "_", OUTCOME, "_significant.txt"))
liberal_sign0$clump_tresh <- "0.2"

liberal_sign <- liberal_sign0[liberal_sign0$exposure %in% mr_mrbase0.001_sign$exposure,]

if (plyr::empty(liberal_sign) == TRUE) {

liberal_sign0$tissue <- EXPOSURE_DATA

message2 <- str_c("no genes significant in conservative analysis (r2 = 0.001)")

write.table(liberal_sign0, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_", EXPOSURE_DATA, "_", OUTCOME, "_significant.txt"), sep = "\t", row.names = F)

} else {

liberal_sign <- distinct(rbind(mr_mrbase0.001_sign,liberal_sign))

liberal_sign$tissue <- EXPOSURE_DATA

write.table(liberal_sign, str_c(EXPOSURE_DATA, "_", OUTCOME, "/results/", "full_results_", EXPOSURE_DATA, "_", OUTCOME, "_significant.txt"), sep = "\t", row.names = F)


message2 <- str_c(length(unique(liberal_sign$exposure)), " genes significant in both main (r2 = 0.2) and conservative analysis (r2 = 0.001)")
}

message1 <- str_c(length(unique(liberal_full$exposure)), " genes tested overall")

}

print(message1)
print(message2)
print("mission_complete")
