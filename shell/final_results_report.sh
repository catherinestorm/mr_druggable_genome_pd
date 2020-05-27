while read OUTCOME; do
    export OUTCOME=${OUTCOME}
    echo "exposure,outcome,nsnp,method,beta,se,p,fdr_qval,clump_tresh,tissue" > full_results/full_results_${OUTCOME}.txt

    echo "exposure,outcome,nsnp,method,beta,se,p,fdr_qval,clump_tresh,tissue,druggability_tier,or,or_lci95,or_uci95,beta_lci95,beta_uci95" > full_results/significant_results_${OUTCOME}.txt

    echo "exposure,outcome,nsnp,egger_intercept,egger_intercept_95_ci,egger_intercept_pvalue,cochrans_q,cochrans_q_pval,i2,tissue" > full_results/full_results_${OUTCOME}_qc.txt

    echo "exposure,outcome,nsnp,egger_intercept,egger_intercept_95_ci,egger_intercept_pvalue,cochrans_q,cochrans_q_pval,i2,tissue" > full_results/significant_results_${OUTCOME}_qc.txt

    while read EXPOSURE_DATA; do
        export OUTCOME=${OUTCOME}
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        Rscript final_results_report.R
        wait
    done < exposure_data.txt
    
done < outcomes.txt
