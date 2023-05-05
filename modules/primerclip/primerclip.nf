process PRIMERCLIP {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'

    container 'quay.io/biocontainers/primerclip:0.3.8--h9ee0642_1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'primerclip' }", pattern:"*.sam", mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/stats' : 'primerclip' }", pattern:"*primerclip_runstats.log", mode:'copy'

    input:
      tuple val(sampleID), file(sam)

    output:
      tuple val(sampleID), file("*.sam"), emit: sam
      tuple val(sampleID), file("*primerclip_runstats.log"), emit: log

    script:

    """
    primerclip ${params.masterfile} ${sam} ${sam.baseName}_primerclip.sam
    """
}
