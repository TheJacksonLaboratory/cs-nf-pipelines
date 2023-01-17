process BOWTIE {
    tag "$sampleID"

    cpus 8
    memory {60.GB * task.attempt}
    time {30.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:daafe97'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'bowtie' }", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sampleID), path(fq_read), val(paired_read_num)

    output:
    tuple val(sampleID), file("*.sam"), emit: sam
    tuple val(sampleID), file("*.log"), emit: bowtie_log

    script:
    """
    zcat ${fq_read} | bowtie -p ${task.cpus} -q -a --best --strata --sam -v 3 -m 100 ${params.bowtie_index} - > ${sampleID}_mapped_${paired_read_num}.sam 2> ${sampleID}.bowtie_${paired_read_num}.log 
    """
    // NOTE: This is hard coded to .gz input files. 

    stub:
    """
    touch ${sampleID}.bowtie_${paired_read_num}.log
    touch ${sampleID}_mapped_${paired_read_num}.sam
    """
}


