process DUPHOLD {
    tag "$sampleID"
    
    cpus = 12
    memory 80.GB
    time '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/duphold:0.2.1--h516909a_1'

    input:
    tuple val(sampleID), path(alignment_file), path(alignement_index)
    tuple val(sampleID), path(sv_variants)
    tuple val(sampleID), path(snp_variants)
    tuple val(sampleID), path(snp_variants_index)
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(sampleID), path("*.vcf.gz")   , emit: vcf

    script:

    def snp_annotation = snp_variants ? "--snp ${snp_variants}" : ""
    name="\$(echo \"${sv_variants}\" | sed 's/Sort.vcf//g')"


    if ("${sv_variants}" =~ "delly")
    """
    export DUPHOLD_SAMPLE_NAME=${name}
  
    duphold \\
        --threads ${task.cpus} \\
        --output ${sampleID}_delly.vcf.gz \\
        --vcf ${sv_variants} \\
        --bam ${alignment_file} \\
        --fasta ${fasta} \\
        ${snp_annotation}
    """
    else if ("${sv_variants}" =~ "manta")
    """
    duphold \\
        --threads ${task.cpus} \\
        --output ${sampleID}_manta.vcf.gz \\
        --vcf ${sv_variants} \\
        --bam ${alignment_file} \\
        --fasta ${fasta} \\
        ${snp_annotation}
    """

}
