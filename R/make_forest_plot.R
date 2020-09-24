## this script will put the results and quality control metrics for all significant outcomes into one data frame

# prep
library(dplyr)
library(readr)
library(forestplot)
library(stringr)
library(tidyverse)


format_numbers <- function(number) {
  if (is.numeric(number) == FALSE & is.integer(number) == FALSE) {
    result <- NA
  } else if (number < 0.001) {
    result <- formatC(number, format = "e", digits = 2)
  } else if (number >= 0.001) {
    result <- signif(number, digits = 3)
  }
  return(result)
}




significant_res <- as.data.frame(read_tsv("full_results/significant_genes_results_all_outcomes.txt"))
significant_res <- unique(significant_res)


# make forest plot RISK
discovery <- significant_res[significant_res$outcome == "nalls2014",]
discovery$exposure_tissue <- str_c(discovery$exposure, "_", discovery$tissue)
replication <- significant_res[significant_res$outcome == "nalls2019",]
replication$exposure_tissue <- str_c(replication$exposure, "_", replication$tissue)
discovery <- subset(discovery, discovery$exposure_tissue %in% replication$exposure_tissue)


data <- as.data.frame(subset(discovery, discovery$clump_tresh == "0.2"))
data <- data[order(data$exposure),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
data$or_ci <- str_c(round(data$or, digits = 2), " (", round(data$or_lci95, digits = 2), ", ",round(data$or_uci95, digits = 2), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

forest_data <- data[,c("or", "or_lci95", "or_uci95", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

table_text <- data[, c("exposure","nsnp", "or_ci", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text <- rbind(c("Gene", "No. SNPs", "OR (95% CI)", "Adjusted P"), table_text)


pdf("figures/forest_risk_risk.pdf", width = 30,height = 30, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = .2,
           align = "l",
           mean = forest_data$or,
           lower = forest_data$or_lci95,
           upper = forest_data$or_uci95,
           xlog=TRUE, # log scale
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours
           xlab="PD odds ratio",
           hrzl_lines = list("2" = gpar(lty = 2)),
           txt_gp = fpTxtGp(summary = list(gpar(cex=3, fontface = "plain")),
                            label = gpar(cex=3),
                            ticks = gpar(cex=3),
                            xlab  = gpar(fontface="bold", cex = 3)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=.2, ...)}
)
dev.off()




# make forest plot AGE AT ONSET

  temp <- significant_res[significant_res$outcome == "blauwendraat2019",]

  data <- as.data.frame(subset(temp, temp$clump_tresh == "0.2"))
  data <- data[order(data$exposure),]

  data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
  data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
  data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
  data$beta_ci <- str_c(signif(data$beta, digits = 3), " (", signif(data$beta_lci, digits = 3), ", ",signif(data$beta_uci95, digits = 3), ")")
  data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

  forest_data <- data[,c("beta", "beta_lci95", "beta_uci95", "tissue")]
  forest_data <- forest_data[complete.cases(forest_data),]
  forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

  table_text <- data[, c("exposure","nsnp", "beta_ci", "roundp")]
  #table_text <- table_text[complete.cases(table_text),]
  table_text <- rbind(c("Gene", "No. SNPs", "Beta (95% CI)", "Adjusted P"), table_text)

  pdf("figures/forest_aao.pdf", width = 30,height = 10, onefile=FALSE)
    forestplot(table_text, # columns to include
               graph.pos = 3, # where graph is
               boxsize = .2,
               align = "l",
               mean = forest_data$beta,
               lower = forest_data$beta_lci95,
               upper = forest_data$beta_uci95,
               is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows
               col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours
               xlab="Beta for PD age at onset",
               hrzl_lines = list("2" = gpar(lty = 2)),
               txt_gp = fpTxtGp(summary = list(gpar(cex=3, fontface = "plain")),
                                label = gpar(cex=3),
                                ticks = gpar(cex=3),
                                xlab  = gpar(fontface="bold", cex = 3)),
               fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=.2, ...)}
    )
    dev.off()






# make forest plot PROGRESSION
data_forest_loop <-  significant_res[!grepl("nalls", significant_res$outcome) & !grepl("blauwendraat", significant_res$outcome),]
data <- as.data.frame(subset(data_forest_loop, data_forest_loop$clump_tresh == "0.2"))
#data[data$outcome == "blauwendraat2019","outcome"] <- "Age at onset"
data <- data[order(data$outcome, data$exposure),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
data$beta_ci <- str_c(signif(data$beta, digits = 3), " (", signif(data$beta_lci, digits = 3), ", ",signif(data$beta_uci95, digits = 3), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)
data <- data[order(data$exposure),]

forest_data <- data[,c("beta", "beta_lci95", "beta_uci95", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

table_text <- data[, c("exposure","outcome", "nsnp", "beta_ci", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text$outcome <- sub(".*?_","",  table_text$outcome)
table_text$outcome <- sub("\\_.*","",  table_text$outcome)
table_text <- rbind(c("Gene","Outcome", "No. SNPs", "Beta (95% CI)", "Adjusted P"), table_text)



pdf("figures/forest_all_progression_outcomes.pdf", width = 45,height = 30, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 4, # where graph is
           boxsize = .2,
           align = "l",
           mean = forest_data$beta,
           lower = forest_data$beta_lci95,
           upper = forest_data$beta_uci95,
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours
           xlab="Beta",
           hrzl_lines = list("2" = gpar(lty = 2)),
           txt_gp = fpTxtGp(summary = list(gpar(cex=4, fontface = "plain")),
                            label = gpar(cex=4),
                            ticks = gpar(cex=4),
                            xlab  = gpar(fontface="bold", cex = 4)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=.2, ...)}
)
dev.off()













# PQTL DATA

## risk
data <- read.table("replication_pqtl_nalls2014/results/full_results_liberal_r2_0.2_replication_pqtl_nalls2014.txt", sep = "\t", header = T)

data <- data[order(data$outcome, data$exposure),]
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
data <- TwoSampleMR::generate_odds_ratios(data[,1:7])

data$pqtl_study0 <- gsub(".*\\_", "", data$exposure)
data$pqtl_study[which(data$pqtl_study0 == "emilsson2018")] <- "Emilsson et al. 2018"
data$pqtl_study[which(data$pqtl_study0 == "sun2018")] <- "Sun et al. 2018"
data$pqtl_study[which(data$pqtl_study0 == "suhre2017")] <- "Suhre et al. 2017"

data$gene <- gsub("\\_.*", "", data$exposure)

data$beta_ci <- str_c(signif(data$beta, digits = 3), " (", signif(data$lo_ci, digits = 3), ", ",signif(data$up_ci, digits = 3), ")")
data$or_ci <- str_c(signif(data$or, digits = 3), " (", signif(data$or_lci95, digits = 3), ", ",signif(data$or_uci95, digits = 3), ")")

data$roundp <- lapply(X = data$p, FUN = format_numbers)
data <- data[order(data$exposure),]

forest_data <- data[,c("or", "or_lci95", "or_uci95")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA),forest_data)

table_text <- data[, c("gene","pqtl_study", "nsnp", "or_ci", "roundp")]
table_text <- rbind(c("Gene","pQTL Source", "No. SNPs", "OR (95% CI)", "Unadjusted P"), table_text)




pdf("figures/forest_pqtl_risk.pdf", width = 40,height = 20, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 4, # where graph is
           boxsize = .2,
           align = "l",
           mean = forest_data$or,
           lower = forest_data$or_lci95,
           upper = forest_data$or_uci95,
           xlog=TRUE,
           col=fpColors(box="darkred", line = "darkred"), # colours
           xlab="PD odds ratio",
           hrzl_lines = list("2" = gpar(lty = 2)),
           txt_gp = fpTxtGp(label = gpar(cex=4),
                            ticks = gpar(cex=4),
                            xlab  = gpar(fontface="bold", cex = 4))
)
dev.off()









## updrs part 4

data <- read.table("replication_pqtl_cont_UPDRS4_scaled/results/full_results_liberal_r2_0.2_replication_pqtl_cont_UPDRS4_scaled.txt", sep = "\t", header = T)


data <- data[order(data$outcome, data$exposure),]
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
data <- TwoSampleMR::generate_odds_ratios(data[,1:7])

data$pqtl_study0 <- gsub(".*\\_", "", data$exposure)
data$pqtl_study[which(data$pqtl_study0 == "emilsson2018")] <- "Emilsson et al. 2018"
data$pqtl_study[which(data$pqtl_study0 == "sun2018")] <- "Sun et al. 2018"
data$pqtl_study[which(data$pqtl_study0 == "suhre2017")] <- "Suhre et al. 2017"

data$gene <- gsub("\\_.*", "", data$exposure)

data$beta_ci <- str_c(signif(data$beta, digits = 3), " (", signif(data$lo_ci, digits = 3), ", ",signif(data$up_ci, digits = 3), ")")
data$or_ci <- str_c(signif(data$or, digits = 3), " (", signif(data$or_lci95, digits = 3), ", ",signif(data$or_uci95, digits = 3), ")")

data$roundp <- lapply(X = data$p, FUN = format_numbers)
data <- data[order(data$exposure),]


forest_data <- data[,c("beta", "lo_ci", "up_ci")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA),forest_data)

table_text <- data[, c("gene","pqtl_study", "nsnp", "beta_ci", "roundp")]
table_text <- rbind(c("Gene","pQTL Source", "No. SNPs", "Beta (95% CI)", "Unadjusted P"), table_text)



pdf("figures/forest_pqtl_updrs.pdf", width = 40,height = 10, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 4, # where graph is
           boxsize = .2,
           align = "l",
           mean = forest_data$beta,
           lower = forest_data$lo_ci,
           upper = forest_data$up_ci,
           col=fpColors(box="darkred", line = "darkred"), # colours
           xlab="UPDRS part 4 score",
           hrzl_lines = list("2" = gpar(lty = 2)),
           txt_gp = fpTxtGp(label = gpar(cex=4),
                            ticks = gpar(cex=4),
                            xlab  = gpar(fontface="bold", cex = 4))
)
dev.off()
