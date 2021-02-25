
library("coloc")
library("stringr")

#GENE <- Sys.Getenv("GENE")
#EXPOSURE_DATA <- Sys.Getenv("EXPOSURE_DATA")
#OUTCOME <- Sys.Getenv("OUTCOME")


# set priors
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5

# read in data
dat_harmonized4coloc_files <- list.files(path = "coloc", pattern = "dat_harmonized4coloc_pd"
,full.names=T)


coloc_results_gene <- list()

for (i in 1:length(dat_harmonized4coloc_files)){

    dat <- read.table(dat_harmonized4coloc_files[i],header=T,sep="\t",stringsAsFactors=F)

    GENE <- dat$exposure[1]
    OUTCOME <- dat$outcome[1]
    EXPOSURE_DATA <- dat$tissue[1]

    dat$var.exposure <- dat$se.exposure^2
    dat$var.outcome <- dat$se.outcome^2


    if (OUTCOME == "pd_risk_discovery") {

        PD_cc_ratio <- 13708/(13708+95282)

        # using eaf as maf

        coloc_results_prior <-
      coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                 varbeta = dat$var.exposure,  # standard error squared
                                 N = dat$samplesize.exposure,  # sample size # might work without it
                                 type = "quant",
                                 sdY = 1,
                                 #MAF = dat$eaf.exposure,
                                 snp = dat$SNP),
                dataset2 = list( beta = dat$beta.outcome,
                                 varbeta = dat$var.outcome,
                                 type = "cc", # case control
                                 N = dat$samplesize.outcome,
                                 s = PD_cc_ratio,
                                 MAF = dat$eaf.outcome,
                                 snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
                             )
    } else if (OUTCOME == "pd_risk_replication") {

        PD_cc_ratio <- 8036/(8036+5803)

        # using eaf as maf
        coloc_results_prior <-
        coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                 varbeta = dat$var.exposure,  # standard error squared
                                 N = dat$samplesize.exposure,  # sample size # might work without it
                                 type = "quant",
                                 sdY = 1,
                                 #MAF = dat$eaf.exposure,
                                 snp = dat$SNP),
                dataset2 = list( beta = dat$beta.outcome,
                                 varbeta = dat$var.outcome,
                                 type = "cc", # case control
                                 N = dat$samplesize.outcome,
                                 s = PD_cc_ratio,
                                 MAF = dat$eaf.outcome,
                                 snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
                             )

    } else if (OUTCOME == "pd_age_at_onset") {

    # using eaf as maf
    coloc_results_prior <-
      coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                 varbeta = dat$var.exposure,  # standard error squared
                                 N = dat$samplesize.exposure,  # sample size # might work without it
                                 type = "quant",
                                 sdY = 1,
                                 #MAF = dat$eaf.exposure,
                                 snp = dat$SNP),
                dataset2 = list( beta = dat$beta.outcome,
                                 varbeta = dat$var.outcome,
                                 N = dat$samplesize.outcome,
                                 type = "quant",
                                 MAF = dat$eaf.outcome,
                                 snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
                             )
    } else if (grepl("cont",OUTCOME) | grepl("surv",OUTCOME)) {

        coloc_results_prior <-
          coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                     varbeta = dat$var.exposure,  # standard error squared
                                     N = dat$samplesize.exposure,  # sample size # might work without it
                                     type = "quant",
                                     sdY = 1,
                                     #MAF = dat$eaf.exposure,
                                     snp = dat$SNP),
                    dataset2 = list( beta = dat$beta.outcome,
                                     varbeta = dat$var.outcome,
                                     N = dat$samplesize.outcome,
                                     type = "quant",
                                     MAF = dat$eaf.outcome,
                                     snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
                                 )
        }


    # make a pretty table

    coloc_results_prior[["summary"]]["p1"] <- coloc_results_prior[["priors"]]["p1"]
    coloc_results_prior[["summary"]]["p2"] <- coloc_results_prior[["priors"]]["p2"]
    coloc_results_prior[["summary"]]["p12"] <- coloc_results_prior[["priors"]]["p12"]

    coloc_results_prior[["summary"]]["gene"] <- GENE
    coloc_results_prior[["summary"]]["outcome"] <- OUTCOME
    coloc_results_prior[["summary"]]["exposure_data"] <- EXPOSURE_DATA

    if(i == 1){

      coloc_results_gene_summary <- coloc_results_prior[["summary"]] %>% t() %>% tibble::as_tibble()

    } else {

      coloc_results_gene_summary <-
        dplyr::bind_rows(coloc_results_gene_summary,
                  coloc_results_prior[["summary"]] %>% t() %>% tibble::as_tibble())

    }
}

# keep only powered analyses
coloc_results_gene_summary$sum_pph3_pph4 <- as.numeric(coloc_results_gene_summary$PP.H3.abf) + as.numeric(coloc_results_gene_summary$PP.H4.abf)

coloc_results_gene_summary_powered <- coloc_results_gene_summary[which(coloc_results_gene_summary$sum_pph3_pph4 > 0.8),]

write.table(coloc_results_gene_summary_powered, str_c("coloc/res_coloc_all_outcomes_p1_",p1,"_p2_",p2,"_p12_",p12, ".txt"),sep="\t",quote=F,row.names=F)

print("mission_complete")
