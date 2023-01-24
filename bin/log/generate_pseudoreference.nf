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
--snp_vcf                       ${params.snp_vcf}
--indel_vcf                     ${params.indel_vcf}
--primary_reference_fasta       ${params.primary_reference_fasta}
--primary_reference_gtf         ${params.primary_reference_gtf}
--strain                        ${params.strain}
--genome_version                ${params.genome_version}
--diploid                       ${params.diploid}
--keep_fails                    ${params.keep_fails}
--pass_only                     ${params.pass_only}
--quality_filter                ${params.quality_filter}
--region                        ${params.region}
--bed                           ${params.bed}

Project Directory: ${projectDir}
______________________________________________________
"""
}