process JAFFA {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 10.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'quay.io/jaxcompsci/jaffa:d1587c9'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/fusions': 'jaffa' }", pattern: "*jaffa_results.csv", mode:'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/fusions': 'jaffa' }", pattern: "*jaffa_results.fasta", mode:'copy', enabled: params.keep_intermediate

    input:
        tuple val(sampleID), path(reads)

    output:
        tuple val(sampleID), path("*jaffa_results.csv"), emit: jaffa_fusions
        tuple val(sampleID), path("*jaffa_results.fasta"), emit: jaffa_fasta

    script:
    ext = reads[0].getExtension()

    """

    bpipe run -v \
    -n ${task.cpus} \
    -p fastqInputFormat='%_*.${ext}' \
    -p refBase=${params.jaffa_ref_dir} \
    -p genome=hg38 \
    -p annotation=genCode22 \
    /opt/JAFFA/JAFFA_direct.groovy \
    ${reads[0]} \
    ${reads[1]}

    mv jaffa_results.csv ${sampleID}_jaffa_results.csv
    mv jaffa_results.fasta ${sampleID}_jaffa_results.fasta ; 

    """
}