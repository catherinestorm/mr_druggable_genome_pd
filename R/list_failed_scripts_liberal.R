library("stringr")

all <- read.table("all_liberal_scripts.txt", colClasses = "character", header = F)
all$outcome <- gsub(".*liberal_r2_0.2_", "", all$V1)
all$outcome <- gsub(".sh", "", all$outcome)


complete <- read.table("completed_liberal_scripts.txt", colClasses = "character", header = F)
complete$outcome <- gsub(".*liberal_r2_0.2_egger_cochransq_", "", complete$V1)
complete$outcome <- gsub(".txt", "", complete$outcome)

not_completed <- all[!(all$outcome %in%complete$outcome),]
not_completed$run <- gsub(".sh", "", not_completed$V1)


write.table(not_completed$run, "failed_scripts_liberal.txt", col.names=F, row.names=F, quote=F)

