process VCFTOOLS_FILTER {
    tag "$sampleID"

    cpus = 1
    memory = 10.GB
    time = '23:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    input:
        tuple val(sampleID), file(merged_vcf)

    output:
        tuple val(sampleID), file("${sampleID}_mergedCall.NS.main.vcf"), emit: vcf

    container 'quay.io/jaxcompsci/vcftools:0.1.17--g581c231'

    script:
        """
        vcftools \
            --vcf ${merged_vcf} \
            --chr 1 \
            --chr 2 \
            --chr 3 \
            --chr 4 \
            --chr 5 \
            --chr 6 \
            --chr 7 \
            --chr 8 \
            --chr 9 \
            --chr 10 \
            --chr 11 \
            --chr 12 \
            --chr 13 \
            --chr 14 \
            --chr 15 \
            --chr 16 \
            --chr 17 \
            --chr 18 \
            --chr 19 \
            --chr MT \
            --chr X \
            --chr Y \
            --chr chr1 \
            --chr chr2 \
            --chr chr3 \
            --chr chr4 \
            --chr chr5 \
            --chr chr6 \
            --chr chr7 \
            --chr chr8 \
            --chr chr9 \
            --chr chr10 \
            --chr chr11 \
            --chr chr12 \
            --chr chr13 \
            --chr chr14 \
            --chr chr15 \
            --chr chr16 \
            --chr chr17 \
            --chr chr18 \
            --chr chr19 \
            --chr chrM \
            --chr chrX \
            --chr chrY \
            --recode \
            --recode-INFO-all \
            --out ${sampleID}_mergedCall.NS.main;
        grep -v "##contig=<ID=GL" ${sampleID}_mergedCall.NS.main.recode.vcf | \
        grep -v "##contig=<ID=JH" > ${sampleID}_mergedCall.NS.main.vcf;
        """
}