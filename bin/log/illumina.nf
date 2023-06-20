def PARAM_LOG(){
    log.info """
______________________________________________________

              ILLUMINA SV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--genome_build         ${params.genome_build}
--fasta                ${params.fasta}
--genome               ${params.genome}
--bwa_index            ${params.bwa_index}
--sample_folder        ${params.sample_folder}
--fastq1               ${params.fastq1}
--fastq2               ${params.fastq2}
--bam                  ${params.bam}
--names                ${params.names}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--quality_phred        ${params.quality_phred}
--unqualified_perc     ${params.unqualified_perc}
--exclude_regions      ${params.exclude_regions}
--tandem_repeats       ${params.tandem_repeats}
--sv_ins_ref           ${params.sv_ins_ref}
--sv_del_ref           ${params.sv_del_ref}
--sv_inv_ref           ${params.sv_inv_ref}
--reg_ref              ${params.reg_ref}
--genes_bed            ${params.genes_bed}
--exons_bed            ${params.exons_bed}
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