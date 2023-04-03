    process GATK_MERGEMUTECTSTATS {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time '05:00:00'
    // time {3.hour * task.attempt}
    // errorStrategy 'retry' 
    // maxRetries 1

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(list)

    output:
    tuple val(sampleID), file("*.stats"), emit: stats

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    stats = list.collect { "-stats $it" }.join(' ')

    """
    gatk --java-options "-Xmx${my_mem}G" MergeMutectStats \
    ${stats} \
    -O ${sampleID}_merged.stats
    """
    }