def param_log(){
    log.info """
______________________________________________________

                PACBIO PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--fasta                ${params.fasta}
--fastq1               ${params.fastq1}
--names                ${params.names}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--surv_dist            ${params.surv_dist}
--surv_supp            ${params.surv_supp}
--surv_type            ${params.surv_type}
--surv_strand          ${params.surv_strand}
--surv_min             ${params.surv_min}

Project Directory: ${projectDir}
______________________________________________________
"""    
}