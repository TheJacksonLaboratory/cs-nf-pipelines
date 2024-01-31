process SEQUENZA_PILEUP2SEQZ {
    tag "$sampleID"

    cpus 2
    memory 60.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/sequenza:v1'

    input:
    tuple val(sampleID), val(meta), path(normal_pileup), val(normal_name), path(tumor_pileup), val(tumor_name)

    output:
    tuple val(sampleID), val(meta), path("*.clean.seqz"), emit: seqz

    script:
    """
    sequenza-utils bam2seqz -p -gc ${params.sequenza_gc} -n ${normal_pileup} -t ${tumor_pileup} | sequenza-utils seqz_binning -w 50 -s - | gzip > ${sampleID}.seqz.gz

    zcat ${sampleID}.seqz.gz | grep -P 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY|chromosome' | grep -v '_random' > ${sampleID}.clean.seqz
    """
}
