process SNPSIFT_EXTRACT_AND_PARSE {

    // NOTE: This script is for the parsing of the 'SV' pipeline germline annotationed table from snpeff extractfields. 
    //       It is hard coded to the annotations used. 

    tag "$sampleID"

    cpus = 1
    memory = 6.GB
    time = '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.txt", mode:'copy'

    input:
    tuple val(sampleID), file(table)

    output:
    tuple val(sampleID), file("*.txt"), emit: txt

    script:
    """
    python ${projectDir}/bin/pta/split_annotations.py ${table} ${sampleID}_annotated_filtered_final_table.txt
    """
}
