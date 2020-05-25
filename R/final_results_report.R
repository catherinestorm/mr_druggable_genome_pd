library("readr")
suppressMessages(library(dplyr, warn.conflict = FALSE, quietly = TRUE))
library("stringr")


outcomes <- read_tsv("outcomes.txt", col_names=F, col_types = cols())
progression_outcomes <- read_tsv("progression_outcomes.txt", col_names=F, col_types = cols())

druggable_genome <- read_tsv("druggable_genome_new.txt", col_types = cols())

druggable_genome <- druggable_genome[,c("gene_display_label", "priority")]
names(druggable_genome) <- c("exposure", "druggability_tier")

progression_results <- data.frame()
progression_results_qc <- data.frame()
sign_progression_results <- data.frame()
sign_progression_results_qc <- data.frame()

for (i in 1:nrow(outcomes)) {
OUTCOME <- paste(outcomes[i,])

if (OUTCOME %in% progression_outcomes$X1) {


# full data
temp_eqtlgen_progression <- read_tsv(str_c("eqtlgen_",OUTCOME ,"/results/full_results_eqtlgen_",OUTCOME,".txt"), col_types = cols())
temp_psychencode_progression <- read_tsv(str_c("psychencode_",OUTCOME ,"/results/full_results_psychencode_",OUTCOME,".txt"), col_types = cols())
temp_total_progression <- distinct(rbind(temp_eqtlgen_progression, temp_psychencode_progression))
progression_results <- distinct(rbind(progression_results, temp_total_progression))


temp_eqtlgen_qc <- read_tsv(str_c("eqtlgen_",OUTCOME ,"/results/full_results_liberal_r2_0.2_egger_cochransq_eqtlgen_",OUTCOME,".txt"), col_types = cols())
temp_psychencode_qc <- read_tsv(str_c("psychencode_",OUTCOME ,"/results/full_results_liberal_r2_0.2_egger_cochransq_psychencode_",OUTCOME,".txt"), col_types = cols())
temp_total_qc <- distinct(rbind(temp_eqtlgen_qc, temp_psychencode_qc))
progression_results_qc <- distinct(rbind(progression_results_qc, temp_total_qc))




# significant
sign_temp_eqtlgen_progression <- read_tsv(str_c("eqtlgen_",OUTCOME ,"/results/full_results_eqtlgen_",OUTCOME,"_significant.txt"), col_types = cols())
sign_temp_psychencode_progression <- read_tsv(str_c("psychencode_",OUTCOME ,"/results/full_results_psychencode_",OUTCOME,"_significant.txt"), col_types = cols())
sign_temp_total_progression <- distinct(rbind(sign_temp_eqtlgen_progression, sign_temp_psychencode_progression))
sign_progression_results <- distinct(rbind(sign_progression_results, sign_temp_total_progression))

sign_progression_results_qc <- subset(progression_results_qc, progression_results_qc$exposure %in% sign_progression_results$exposure)



} else {

### all results
eqtlgen <- read_tsv(str_c("eqtlgen_",OUTCOME ,"/results/full_results_eqtlgen_",OUTCOME,".txt"), col_types = cols())
psychencode <- read_tsv(str_c("psychencode_",OUTCOME ,"/results/full_results_psychencode_",OUTCOME,".txt"), col_types = cols())
total <- distinct(rbind(eqtlgen, psychencode))

# add druggability info & calculate odds ratio
total <- left_join(total, druggable_genome, by = "exposure")

total$or <- exp(total$beta)
total$or_lci95 <- exp(total$beta-1.96*total$se)
total$or_uci95 <- exp(total$beta+1.96*total$se)

total$beta_lci95 <- total$beta-1.96*total$se
total$beta_uci95 <- total$beta+1.96*total$se

### significant results
sign_eqtlgen <- read_tsv(str_c("eqtlgen_",OUTCOME ,"/results/full_results_eqtlgen_",OUTCOME,"_significant.txt"), col_types = cols())
sign_psychencode <- read_tsv(str_c("psychencode_",OUTCOME ,"/results/full_results_psychencode_",OUTCOME,"_significant.txt"), col_types = cols())
sign_total <- distinct(rbind(sign_eqtlgen, sign_psychencode))

# add druggability info & calculate odds ratio
if (plyr::empty(sign_total) == F) {sign_total <- left_join(sign_total, druggable_genome, by = "exposure")

sign_total$or <- exp(sign_total$beta)
sign_total$or_lci95 <- exp(sign_total$beta-1.96*sign_total$se)
sign_total$or_uci95 <- exp(sign_total$beta+1.96*sign_total$se)

sign_total$beta_lci95 <- sign_total$beta-1.96*sign_total$se
sign_total$beta_uci95 <- sign_total$beta+1.96*sign_total$se
}

# write out results
write.table(total, str_c("full_results/full_results_",OUTCOME,".txt"), sep = "\t", row.names = F)
write.table(sign_total, str_c("full_results/significant_results_",OUTCOME,".txt"), sep = "\t", row.names = F)

# qc
eqtlgen_qc <- read_tsv(str_c("eqtlgen_",OUTCOME ,"/results/full_results_liberal_r2_0.2_egger_cochransq_eqtlgen_",OUTCOME,".txt"), col_types = cols())
psychencode_qc <- read_tsv(str_c("psychencode_",OUTCOME ,"/results/full_results_liberal_r2_0.2_egger_cochransq_psychencode_",OUTCOME,".txt"), col_types = cols())
total_qc <- distinct(rbind(eqtlgen_qc, psychencode_qc))
sign_total_qc <- subset(total_qc, total_qc$exposure %in% sign_total$exposure)

write.table(total_qc, str_c("full_results/full_results_",OUTCOME,"_qc.txt"), sep = "\t", row.names = F)
write.table(sign_total_qc, str_c("full_results/significant_results_",OUTCOME,"_qc.txt"), sep = "\t", row.names = F)

print(str_c(OUTCOME, ": ", length(unique(total$exposure)), " genes tested in total"))
print(str_c(OUTCOME, ": ", length(unique(eqtlgen$exposure)), " genes tested in eqtlgen"))
print(str_c(OUTCOME, ": ", length(unique(psychencode$exposure)), " genes tested in psychencode"))

print("")

print(str_c(OUTCOME, ": ", length(unique(sign_total$exposure)), " genes reaching significance in total"))
print(str_c(OUTCOME, ": ", length(unique(sign_eqtlgen$exposure)), " genes reaching significance in eqtlgen"))
print(str_c(OUTCOME, ": ", length(unique(sign_psychencode$exposure)), " genes reaching significance in psychencode"))

print("")
print("")
}}





print(str_c("progression: ", length(unique(progression_results$exposure)), " genes tested in total"))
print(str_c("progression: ", length(unique(progression_results$exposure[progression_results$tissue == "eqtlgen"])), " genes tested in eqtlgen"))
print(str_c("progression: ", length(unique(progression_results$exposure[progression_results$tissue == "psychencode"])), " genes tested in psychencode"))


print(str_c("progression: ", length(unique(sign_progression_results$exposure)), " genes reaching significance in total"))
print(str_c("progression: ", length(unique(sign_progression_results$exposure[sign_progression_results$tissue == "eqtlgen"])), " genes reaching significance in eqtlgen"))
print(str_c("progression: ", length(unique(sign_progression_results$exposure[sign_progression_results$tissue == "psychencode"])), " genes reaching significance in psychencode"))

print("")



display_progression <- distinct(progression_results[,c("exposure", "outcome")])
display_progression_eqtlgen <- distinct(progression_results[progression_results$tissue == "eqtlgen",c("exposure", "outcome")])
display_progression_psychencode <- distinct(progression_results[progression_results$tissue == "psychencode",c("exposure", "outcome")])

progression <- suppressWarnings(data.frame(table(display_progression$outcome)) %>% 
            left_join(data.frame(table(display_progression_eqtlgen$outcome)), by = "Var1") %>% 
            left_join(data.frame(table(display_progression_psychencode$outcome)), by = "Var1"))
names(progression) <- c("outcome", "tested", "tested_eqtlgen","tested_psychencode")
progression$outcome <- as.character(progression$outcome)
progression <- rbind(progression, c("total", length(unique(display_progression$exposure)), 
                                                length(unique(display_progression_eqtlgen$exposure)), 
                                                length(unique(display_progression_psychencode$exposure))))
progression

print("")


display_progression_significant <- distinct(sign_progression_results[,c("exposure", "outcome")])
display_progression_significant_eqtlgen <- distinct(sign_progression_results[sign_progression_results$tissue == "eqtlgen",c("exposure", "outcome")])
display_progression_significant_psychencode <- distinct(sign_progression_results[sign_progression_results$tissue == "psychencode",c("exposure", "outcome")])

significant <- suppressWarnings(data.frame(table(display_progression_significant$outcome)) %>% 
            left_join(data.frame(table(display_progression_significant_eqtlgen$outcome)), by = "Var1") %>% 
            left_join(data.frame(table(display_progression_significant_psychencode$outcome)), by = "Var1"))
names(significant) <- c("outcome", "significant", "significant_eqtlgen", "significant_psychencode")
significant$outcome <- as.character(significant$outcome)
significant <- rbind(significant, c("total", length(unique(display_progression_significant$exposure)), 
                                                length(unique(display_progression_significant_eqtlgen$exposure)), 
                                                length(unique(display_progression_significant_psychencode$exposure))))
significant



progression_results <- left_join(progression_results, druggable_genome, by = "exposure")
progression_results$or <- exp(progression_results$beta)
progression_results$or_lci95 <- exp(progression_results$beta-1.96*progression_results$se)
progression_results$or_uci95 <- exp(progression_results$beta+1.96*progression_results$se)

progression_results$beta_lci95 <- progression_results$beta-1.96*progression_results$se
progression_results$beta_uci95 <- progression_results$beta+1.96*progression_results$se

write.table(progression_results, str_c("full_results/full_results_progression.txt"), sep = "\t", row.names = F)
write.table(progression_results_qc, str_c("full_results/full_results_progression_qc.txt"), sep = "\t", row.names = F)

sign_progression_results <- left_join(sign_progression_results, druggable_genome, by = "exposure")
sign_progression_results$or <- exp(sign_progression_results$beta)
sign_progression_results$or_lci95 <- exp(sign_progression_results$beta-1.96*sign_progression_results$se)
sign_progression_results$or_uci95 <- exp(sign_progression_results$beta+1.96*sign_progression_results$se)

sign_progression_results$beta_lci95 <- sign_progression_results$beta-1.96*sign_progression_results$se
sign_progression_results$beta_uci95 <- sign_progression_results$beta+1.96*sign_progression_results$se

write.table(sign_progression_results, str_c("full_results/significant_results_progression.txt"), sep = "\t", row.names = F)
write.table(progression_results_qc, str_c("full_results/significant_results_progression_qc.txt"), sep = "\t", row.names = F)
