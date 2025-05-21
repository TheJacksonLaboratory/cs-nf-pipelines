process SEQUENZA_RUN {
    tag "$sampleID"

    cpus 2
    memory 60.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/sequenza:v1'

    publishDir "${params.pubdir}/${sampleID + '/sequenza_cnv/txt'}", pattern:"*.txt", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/sequenza_cnv'}", pattern:"*segments.txt", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/sequenza_cnv/Rdata'}", pattern:"*.RData", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/sequenza_cnv/pdfs'}", pattern:"*.pdf", mode:'copy'

    input:
    tuple val(sampleID), val(meta), path(seqz)

    output:
    tuple val(sampleID), val(meta), path("*extract.RData"), emit: extract_rdata
    tuple val(sampleID), val(meta), path("*table.RData"), emit: table_rdata
    tuple val(sampleID), val(meta), path("*segments.txt"), emit: segments
    tuple val(sampleID), path("*_segments.tmp"), emit: segments_tmp
    tuple val(sampleID), path("*.pdf"), emit: pdf
    tuple val(sampleID), path("*.txt"), emit: txt
    

    script:    
    female = meta.sex == 'XX' ? 'TRUE' : 'FALSE'

    """
    Rscript ${projectDir}/bin/wes/sequenza_run.R ${seqz} ${sampleID} ./ ${female}

    cat ${sampleID}_segments.txt | tr -d "\\"" | awk 'NR>1' > ${sampleID}_segments.tmp
    """
}
// NOTE: If sample is XX, female, else if XY or NA run sample as male. 
// the `cat` statement prepares the segments file for bedtools substraction against the NA windows file. 

// Note: header of the Sequenza segments file is used in the bedtools substract step. If header is ever altered by sequenza devs, the header file will require an update.
