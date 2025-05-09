process DUPHOLD {
    tag "$sampleID"
    
    cpus = 12
    memory 80.GB
    time '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/duphold:0.2.1--h516909a_1'

    publishDir "${params.pubdir}/${sampleID + '/duphold'}", pattern:"*.vcf", mode:'copy'

    input:
    tuple val(sampleID), path(alignment_file), path(alignement_index), path(sv_variants)
    tuple path(fasta), path(fasta_fai)
    val (tool)

    output:
    tuple val(sampleID), path("*.vcf"), emit: vcf

    script:
    """
    export DUPHOLD_SAMPLE_NAME=${sampleID}_${tool}
  
    duphold \
    --threads ${task.cpus} \
    --output ${sampleID}_${tool}_duphold.vcf.gz \
    --vcf ${sv_variants} \
    --bam ${alignment_file} \
    --fasta ${fasta}
   
    gunzip ${sampleID}_${tool}_duphold.vcf.gz

    """
}
