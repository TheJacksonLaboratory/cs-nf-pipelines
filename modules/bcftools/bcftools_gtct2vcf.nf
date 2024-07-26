def prepare_bcftools_inputs(ch_input) {
    bpm_file = file(params.bpm_file)
    csv_file = file(params.csv_input)
    egt_file = file(params.egt_file)
    gtcs_dir = IAAP_CLI.out.gtc
    fasta_file = file(params.fasta_file)
    tsv_file = file(params.tsv_file)

    return tuple(bpm_file, csv_file, egt_file, gtcs_dir, fasta_file, tsv_file)
}


// Define BCFTOOLS_GTC2VCF process
process BCFTOOLS_GTC2VCF {
    cpus = 4
    memory 24.GB
    time '01:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'
    publishDir "${params.pubdir}", mode: 'copy'

    input:
    tuple path(bpm_file), path(csv_file), path(egt_file), path(gtcs_dir), path(fasta_file), path(tsv_file) from prepare_bcftools_inputs

    output:
    tuple path('bcftools_convert.bcf'), path('bcftools_convert.bcf.csi'), path('bcftools_convert.vcf'), path('bcftools_convert.tsv')

    script:
    """
    bcftools +gtc2vcf --no-version -Ou \
    --bpm ${bpm} \
    --csv ${csv} \
    --egt ${egt} \
    --gtcs ${gtcs_dir} \
    --fasta-ref ${fasta} \
    --extra ${tsv} | \
    bcftools sort -Ou -T ./bcftools. | \
    bcftools norm --no-version -Ob -c x -f ${fasta} | \
    tee bcftools_convert.bcf | \
    bcftools index --force --output bcftools_convert.bcf.csi
    bcftools convert -O v -o bcftools_convert.vcf bcftools_convert.bcf
    """
}
