process IAAP_CLI {
    tag "$sampleID"
    
    cpus = 1
    memory 24.GB
    time '01:30:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.gtc", mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.log", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), path(red_idat), path(green_idat)

    output:
    tuple val(sampleID), val(meta), path("*.gtc"), emit: gtc
    tuple val(sampleID), val(meta), emit: ascat2r
    path "iaap_cli.log", emit: iaap_cli_log

    script:
    """
    /usr/local/bin/iaap-cli/iaap-cli gencall \
        ${params.bpm_file} \
        ${params.egt_file} \
        ./ \
        --idat-folder ./ \
        --output-gtc >> iaap_cli.log 2>&1
    """
}
