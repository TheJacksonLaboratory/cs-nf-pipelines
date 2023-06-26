process FILTER_GTF {
    tag "Filtering GTF"

    cpus 1
    memory 5.GB
    time 2.hour

    container 'quay.io/jaxcompsci/ripgrep:v1'

    publishDir "${params.pubdir}", pattern: '*.gtf', mode:'copy'

    output:
    path('*.gtf'), emit: filtered_gtf

    script:
    include_statement = params.gtf_biotype_include.split(',').collect { "$it" }.join('|')
    """
    sh ${projectDir}/bin/g2gtools/filter_gtf.sh ${params.primary_reference_gtf} "${include_statement}"
    """

    stub:
    """
    touch "gtf.included.gtf"
    """

}
