def PARAM_LOG(){
    log.info """
______________________________________________________

              PACBIO SV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--genome_build         ${params.genome_build}
--pbmode               ${params.pbmode}
--fasta                ${params.fasta}
--fastq1               ${params.fastq1}
--bam                  ${params.bam}
--sampleID             ${params.sampleID}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--pbsv_tandem          ${params.pbsv_tandem}
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