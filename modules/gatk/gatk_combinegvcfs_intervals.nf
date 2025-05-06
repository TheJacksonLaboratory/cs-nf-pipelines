process GATK_COMBINEGVCFS_INTERVALS {
    tag "$sampleID"
    
    cpus 1
    array 23
    memory 16.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    // publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.gvcf", mode:'copy'

    input:
    tuple val(sampleID), path(gvcf), path(idx), val(chrom)

    output:
    tuple val(sampleID), path("*.gvcf"), path("*.idx"), emit: gvcf_idx

    script:
    // memory needs to be set explicitly
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = gvcf.collect { "--variant $it" }.join(' ')

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" CombineGVCFs \
    -R ${params.ref_fa} \
    ${inputs} \
    --intervals ${chrom} \
    -O ${sampleID}_GATKcombined_${chrom}.gvcf
    """
}
