process MULTIQC {

    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/multiqc' : 'multiqc' }", pattern: "*multiqc_report.html", mode:'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/multiqc' : 'multiqc' }", pattern: "*_data", mode:'copy'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data" , emit: data
    path "*_plots" , optional:true, emit: plots

    script:

    """
    multiqc -f $args .
    """

}