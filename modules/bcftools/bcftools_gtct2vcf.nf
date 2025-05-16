process BCFTOOLS_GTC2VCF {
    tag "$sampleID"
    
    cpus = 1
    memory 24.GB
    time '01:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.{vcf,tsv}", mode: 'copy'

    input:
    tuple val(sampleID), val(meta), path(gtc)

    output:
    tuple val(sampleID), val(meta), path('*_convert.bcf'), path('*_convert.bcf.csi'), path('*_convert.vcf'), path('*_convert_info.tsv'), emit: gtc2vcf

    script:
    """
    bcftools +gtc2vcf --no-version -Ou \
    --bpm ${params.bpm_file} \
    --csv ${params.gtc_csv} \
    --egt ${params.egt_file} \
    --gtcs ./ \
    --fasta-ref ${params.ref_fa} \
    --extra ${sampleID}_convert_info.tsv | \
    bcftools sort -Ou -T ./bcftools. | \
    bcftools norm --no-version -Ob -c x -f ${params.ref_fa} | \
    tee ${sampleID}_convert.bcf | \
    bcftools index --force --output ${sampleID}_convert.bcf.csi
    bcftools convert -O v -o ${sampleID}_convert.vcf ${sampleID}_convert.bcf
    """
}
