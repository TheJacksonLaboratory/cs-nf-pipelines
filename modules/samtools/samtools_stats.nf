process SAMTOOLS_STATS {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/samtools' : 'samtools'}"
    }, pattern: "*.flagstat", mode: 'copy'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/samtools' : 'samtools'}"
    }, pattern: "*.idxstats", mode: 'copy'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/samtools' : 'samtools'}"
    }, pattern: "*.stat", mode: 'copy'


    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.flagstat"), emit: flagstat
    tuple val(sampleID), file("*.idxstats"), emit: idxstat
    tuple val(sampleID), file("*.stats"), emit: stats

    script:

    """
    samtools flagstat ${bam[0]} > ${bam[0]}.flagstat
    samtools idxstats ${bam[0]} > ${bam[0]}.idxstats
    samtools stats ${bam[0]} > ${bam[0]}.stats
    """
}
