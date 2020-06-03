

library(readr)




# functions

## for qqplot, adapted from qqman package
my_qq <- function(pvector, ...) {
  
  o = -log10(sort(as.numeric(pvector)))
  e = -log10(ppoints(length(as.numeric(pvector))))
  
  ##define  default arguments
  def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)), 
                   xlab=expression(Expected~~-log[10](italic(p))), 
                   ylab=expression(Observed~~-log[10](italic(p)))
  )
  ##get a list of "..." arguments
  dotargs <- list(...)
  
  ##call the plot function passing NA, your ... arguments, and the default
  tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)
  
  #add diagonal
  abline(0,1,col="blue")
}





# read in meta file
meta <- read_tsv("outcome_data/METAANALYSIS_samplesize_1.tbl")



# quality control

## QQ plot
jpeg('outcome_data/qq_meta_samplesize1.jpeg')
qq(meta$`P-value`, main = "Q-Q Plot: GWAS meta-analysis")
dev.off()



## LAMBDA
chisq <- qchisq(meta$`P-value`,1)
rawLambda <- median(chisq)/qchisq(0.5,1)
lambda1000 <- 1+(rawLambda-1)*(((1/8036)+(1/5803))/((1/1000)+(1/1000)))

lambda_summary <- data.frame("rawLambda" = rawLambda, "lambda1000" = lambda1000)

write.table(lambda_summary, "outcome_data/meta_lambda_summary.txt", row.names = T, sep = "\t")

