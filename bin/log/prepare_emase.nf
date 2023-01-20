def param_log(){
log.info """
______________________________________________________

    GENERATE MULTIWAY TRANSCRIPTOME PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--genome_file_list              ${params.genome_file_list}
--gtf_file_list                 ${params.gtf_file_list}
--haplotype_list                ${params.haplotype_list}

Project Directory: ${projectDir}
______________________________________________________
"""
}