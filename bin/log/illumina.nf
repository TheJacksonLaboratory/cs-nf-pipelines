def PARAM_LOG(){
    log.info """
______________________________________________________

              ILLUMINA SV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--fasta                ${params.fasta}
--fastq1               ${params.fastq1}
--fastq2               ${params.fastq2}
--bam                  ${params.bam}
--names                ${params.names}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--exclude_regions      ${params.exclude_regions}
--quality_phred        ${params.quality_phred} 
--unqualified_perc     ${params.unqualified_perc}
--surv_dist            ${params.surv_dist}
--surv_supp            ${params.surv_supp}
--surv_type            ${params.surv_type}
--surv_strand          ${params.surv_strand}
--surv_min             ${params.surv_min}
--keep_intermediate    ${params.keep_intermediate}

Project Directory: ${projectDir}
______________________________________________________
"""
}