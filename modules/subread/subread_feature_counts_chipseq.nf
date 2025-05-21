process SUBREAD_FEATURECOUNTS {
    tag "${antibody}"

    cpus 4
    memory 4.GB
    time '10:00:00'
  
    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'

    publishDir "${params.pubdir}/${'consensusCalling_' + antibody + '/subread'}", pattern: "*.txt*", mode: 'copy'  

    input:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), path(bams), path(saf)

    output:
    tuple val(antibody), file("*featureCounts.txt")        , emit: counts
    tuple val(antibody), file("*featureCounts.txt.summary"), emit: summary

    script:
    prefix = "${antibody}.consensus_peaks"
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
    pe_params = params.read_type == 'SE' ? '' : '-p --donotsort'
    """
    featureCounts \\
        -F SAF \\
        -O \\
        --fracOverlap 0.2 \\
        -T $task.cpus \\
        $pe_params \\
        -a $saf \\
        -o ${prefix}.featureCounts.txt \\
        ${bam_files.join(' ')}
    """
}
