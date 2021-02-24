
library(readr)
library(stringr)
library(TwoSampleMR)


# signficant genes from MR to test with coloc

genes4coloc_full <- read_tsv("coloc/significant_genes_results_all_outcomes.txt")

genes4coloc <- dplyr::distinct(genes4coloc_full[,c("exposure","outcome","tissue")])


########### EXPOSURES ###########


# eqtlgen

eqtlgen_full <- read_exposure_data(
  filename = "eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  phenotype_col = "GeneSymbol",
  samplesize_col = "NrSamples",
  min_pval = 1e-400
)


## keep significant genes only

eqtlgen <- subset(eqtlgen_full, eqtlgen_full$exposure %in% genes4coloc$exposure[which(genes4coloc$tissue=="eqtlgen")])



# psychencode

psychencode_full <- read_exposure_data(
  filename = "eqtl_data_psychencode/psychencode_exposure_dat_snps_5kb_window.txt",
  sep = "\t",
  snp_col = "Rsid",
  beta_col = "regression_slope",
  se_col = "se",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "nominal_pval",
  phenotype_col = "gene.exposure",
  min_pval = 1e-400
)

## add sample size

psychencode_full$samplesize.exposure <- 1387

## keep significant genes only

psychencode <- subset(psychencode_full, psychencode_full$exposure %in% genes4coloc$exposure[which(genes4coloc$tissue=="psychencode")])


########### OUTCOMES ###########




# risk
for(EXPOSURE_DATA in c("eqtlgen","psychencode")) {

    for(OUTCOME in unique(genes4coloc$outcome)) {

        name <- str_c("dat_harmonized_",OUTCOME, "_", EXPOSURE_DATA)

        genes_to_replicate <- subset(genes4coloc, genes4coloc$outcome==OUTCOME & genes4coloc$tissue==EXPOSURE_DATA)

        exp <- eval(parse(text = EXPOSURE_DATA))

        exp <- subset(exp, exp$exposure %in% genes_to_replicate$exposure)


        if (plyr::empty(exp)==FALSE) {



        if (OUTCOME == "pd_risk_discovery") {

            out <- read_outcome_data(snps = exp$SNP,
                                     filename = "outcome_data/pd_risk_discovery_discovery_risk.txt",
                                     sep = "\t",
                                     snp_col = str_c("rsid_", EXPOSURE_DATA),
                                     beta_col = "Effect",
                                     se_col = "StdErr",
                                     effect_allele_col = "Allele1",
                                     other_allele_col = "Allele2",
                                     eaf_col = "Freq1",
                                     pval_col = "P.value")
            out$outcome <- "pd_risk_discovery"

            out$samplesize.outcome <- 108990

        } else if (OUTCOME == "pd_risk_replication") {

            out <- read_outcome_data(snps = exp$SNP,
                                                 filename = "outcome_data/pd_replication_replication_risk.txt",
                                                 sep = "\t",
                                                 snp_col = str_c("rsid_", EXPOSURE_DATA),
                                                 beta_col = "beta",
                                                 se_col = "se",
                                                 effect_allele_col = "Allele1",
                                                 other_allele_col = "Allele2",
                                                 eaf_col = "Freq1",
                                                 pval_col = "P-value",
                                                 samplesize_col = "Weight")

            out$outcome <- "pd_risk_replication"


        } else if (OUTCOME == "pd_age_at_onset") {

            out <- read_outcome_data(filename = "outcome_data/pd_age_at_onset_replication_aao.txt",
                                             sep = "\t",
                                             snp_col = str_c("rsid_", EXPOSURE_DATA),
                                             beta_col = "Effect",
                                             se_col = "StdErr",
                                             effect_allele_col = "Allele1",
                                             other_allele_col = "Allele2",
                                             eaf_col = "Freq1",
                                             pval_col = "P-value")

            out$outcome <- "pd_age_at_onset"

            out$samplesize.outcome <- 17996


        } else if (grepl("pd_progression",OUTCOME)) {

            OUTCOME_SHORT <- gsub("pd_progression_","",OUTCOME)

            out <- read_outcome_data(
                             filename = str_c("outcome_data/",OUTCOME_SHORT,"_with_alleles.txt"),
                             sep = "\t",
                             snp_col = str_c("rsid_", EXPOSURE_DATA),
                             beta_col = "BETA",
                             se_col = "SE",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "P",
                             eaf_col = "maf",
                             samplesize_col = "N",
                             phenotype_col = "outcome"
                             )
        }


        harmonized <- harmonise_data(
            exposure_dat = exp,
            outcome_dat = out)

        harmonized$tissue <- EXPOSURE_DATA


        assign(name, harmonized)

        for (GENE in unique(harmonized$exposure)) {

            temp <- subset(harmonized, harmonized$exposure==GENE)

            write.table(temp, str_c("coloc/dat_harmonized4coloc_",OUTCOME,"_",EXPOSURE_DATA,"_",GENE,".txt"),sep="\t",quote=F,row.names=F)
        }

}}}



# put all into one data frame

dat4rbind <- ls(pat="dat_harmonized_")

full_data <- data.frame()

for(DAT in dat4rbind) {
    temp <- eval(parse(text = DAT))
    temp <- temp[,!(names(temp) == "chr.exposure")]
    full_data <- rbind(full_data, temp)
}


full_data$outcome <- gsub("pd_progression_cont_","pd_progression_pd_progression_cont_",full_data$outcome)
full_data$outcome <- gsub("pd_progression_surv_","pd_progression_pd_progression_surv_",full_data$outcome)


# sanity check
full_data$exposure_outcome_tissue <- str_c(full_data$exposure,"_",full_data$outcome,"_",full_data$tissue)

genes4coloc$exposure_outcome_tissue <- str_c(genes4coloc$exposure,"_",genes4coloc$outcome,"_",genes4coloc$tissue)


all(sort(unique(genes4coloc$exposure_outcome_tissue)) == sort(unique(full_data$exposure_outcome_tissue))
)

write.table(full_data, "coloc/dat_harmonized4coloc_full.txt",sep="\t",quote=F,row.names=F)

print("mission_complete")
