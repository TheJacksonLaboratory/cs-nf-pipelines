process FUSION_REPORT {
    
    tag "$sampleID"

    cpus 1
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy 'finish'

    container 'quay.io/biocontainers/fusion-report:2.1.5--pyhdfd78af_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/Fusion-Report/' : 'star-fusion' }", mode:'copy'

    input:
        //set val(sample), file(reads), file(arriba), file(ericscript), file(fusioncatcher), file(pizzly), file(squid), file(starfusion) from files_and_reports_summary
        // code left if ever expanded to multiple tools. 

        tuple val(sampleID), file(starfusion)

    output:
        tuple val(sampleID), file("${sampleID}_fusion_list.tsv"), emit: fusion_inspector_input_list
        tuple val(sampleID), file("${sampleID}_fusion_genes_mqc.json"), emit: summary_fusions_mq
        tuple val(sampleID), file("*"), emit: report
    
    script:
    def extra_params = params.fusion_report_opt ? params.fusion_report_opt : ''
    // def tools = !arriba.empty() ? "--arriba ${arriba} " : ''
    // tools += !ericscript.empty() ? "--ericscript ${ericscript} " : ''
    // tools += !fusioncatcher.empty() ? "--fusioncatcher ${fusioncatcher} " : ''
    // tools += !pizzly.empty() ? "--pizzly ${pizzly} " : ''
    // tools += !squid.empty() ? "--squid ${squid} " : ''
    // tools += !starfusion.empty() ? "--starfusion ${starfusion} " : ''
    // Code left if ever expanded to multiple tools. 

    """
    fusion_report run ${sampleID} . ${params.databases} --starfusion ${starfusion} ${extra_params}
    mv fusion_list.tsv ${sampleID}_fusion_list.tsv
    mv fusion_list_filtered.tsv ${sampleID}_fusion_list_filtered.tsv
    mv fusion_genes_mqc.json ${sampleID}_fusion_genes_mqc.json
    """
}