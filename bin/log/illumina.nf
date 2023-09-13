def PARAM_LOG(){
    if (!params.fasta || params.fasta == "<PATH>" || params.fasta == "/<PATH>") {
        error "'--fasta': \"${params.fasta}\" is not valid, specify path to reference FASTA" 
    }

    if (!params.sampleID || params.fasta == "<STRING>") {
        error "'--sampleID': \"${params.sampleID}\" is not valid, specify a name for this sample" 
    }
    
    log.info """
______________________________________________________

              ILLUMINA SV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--data_type            ${params.data_type}
--genome_build         ${params.genome_build}
--fasta                ${params.fasta}
--bwa_index            ${params.bwa_index}
--sample_folder        ${params.sample_folder}
--fastq1               ${params.fastq1}
--fastq2               ${params.fastq2}
--bam                  ${params.bam}
--sampleID             ${params.sampleID}
--pubdir               ${params.pubdir}
-w                     ${workDir}
-c                     ${params.config}
--quality_phred        ${params.quality_phred}
--unqualified_perc     ${params.unqualified_perc}
--exclude_regions      ${params.exclude_regions}
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
