process SURVIVOR_MERGE {

    tag "$sampleID"

    cpus 8
    memory 40.GB
    time "4:00:00"

    container 'quay.io/biocontainers/survivor:1.0.7--he513fc3_0'

    input:
        tuple val(sampleID), file(vcf_tuple)
    output:
        tuple val(sampleID), file("${sampleID}_mergedCall.*.vcf"), emit: vcf
    script:
        if (params.workflow == "pacbio")
            """
            ls ${vcf_tuple[0]} > vcf_list.txt
            ls ${vcf_tuple[1]} >> vcf_list.txt
            SURVIVOR merge vcf_list.txt ${params.surv_dist} ${params.surv_supp} ${params.surv_type} ${params.surv_strand} 0 ${params.surv_min} ${sampleID}_mergedCall.PS.vcf
            """
        else if (params.workflow == "illumina")
            """
            ls ${vcf_tuple[0]} > vcf_list.txt
            ls ${vcf_tuple[1]} >> vcf_list.txt
            ls ${vcf_tuple[2]} >> vcf_list.txt
            ls ${vcf_tuple[3]} >> vcf_list.txt
            SURVIVOR merge vcf_list.txt ${params.surv_dist} ${params.surv_supp} ${params.surv_type} ${params.surv_strand} 0 ${params.surv_min} ${sampleID}_mergedCall.BDLM.vcf
            """
        else if (params.workflow == "ont")
            """
            ls ${vcf_tuple[0]} > vcf_list.txt
            ls ${vcf_tuple[1]} >> vcf_list.txt
            SURVIVOR merge vcf_list.txt ${params.surv_dist} ${params.surv_supp} ${params.surv_type} ${params.surv_strand} 0 ${params.surv_min} ${sampleID}_mergedCall.NS.vcf
            """            
}