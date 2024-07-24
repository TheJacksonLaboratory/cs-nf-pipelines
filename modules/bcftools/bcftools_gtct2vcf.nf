// bcftools_gtc2vcf.nf

process BCFTOOLS_GTC2VCF {
      
    cpus = 4
    memory 24.GB
    time '01:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'
    publishDir "${params.pubdir}", mode: 'copy'

    input:
    tuple path(bpm), path(csv), path(egt), path(gtcs_dir), path(fasta), path(tsv)

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