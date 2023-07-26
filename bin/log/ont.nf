def PARAM_LOG(){
    log.info """
______________________________________________________

            OXFORD NANOPORE SV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--genome_build         ${params.genome_build}
--fasta                ${params.fasta}
--fasta_index          ${params.fasta_index}
--fastq1               ${params.fastq1}
--bam                  ${params.bam}
--sampleID             ${params.sampleID}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--quality              ${params.quality}
--length               ${params.length} 
--headcrop             ${params.headcrop}
--tailcrop             ${params.tailcrop}
--tandem_repeats       ${params.tandem_repeats}
--sv_ins_ref           ${params.sv_ins_ref}
--sv_del_ref           ${params.sv_del_ref}
--sv_inv_ref           ${params.sv_inv_ref}
--reg_ref              ${params.reg_ref}
--genes_bed            ${params.genes_bed}
--exons_bed            ${params.exons_bed}
--targ_chr             ${params.targ_chr}
--targ_start           ${params.targ_start}
--targ_end             ${params.targ_end}
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