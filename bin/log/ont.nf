def PARAM_LOG(){
    log.info """
______________________________________________________

            OXFORD NANOPORE SV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--fasta                ${params.fasta}
--fastq1               ${params.fastq1}
--bam                  ${params.bam}
--names                ${params.names}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--quality              ${params.quality}
--length               ${params.length} 
--headcrop             ${params.headcrop}
--tailcrop             ${params.tailcrop}
--tandem_repeats       ${params.tandem_repeats}
--targ_chr             ${params.targ_chr}
--targ_start           ${params.targ_start}
--targ_end             ${params.targ_end}
--surv_dist            ${params.surv_dist}
--surv_supp            ${params.surv_supp}
--surv_type            ${params.surv_type}
--surv_strand          ${params.surv_strand}
--surv_min             ${params.surv_min}

Project Directory: ${projectDir}
______________________________________________________
"""
}