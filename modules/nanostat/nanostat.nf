process NANOSTAT{
    tag "$sampleID"

    cpus 16
    memory 24.GB
    time "24:00:00"

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'stats'}", pattern: "nanostat*", mode:'copy'

    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'

    input:
        tuple val(sampleID), file(read1)

    output:
        path("nanostat_*"), emit: nanostat_report

    script:
        """
        NanoStat --fastq ${read1} -n nanostat_${read1.baseName}_${sampleID}
        """

}