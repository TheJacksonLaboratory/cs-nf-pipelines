process GATK_VARIANTFILTRATION {
    tag "$sampleID"

    cpus = 1
    memory = 6.GB
    time = '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), file(vcf), file(idx)
    val(output_suffix)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf
    tuple val(sampleID), file("*.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" VariantFiltration \
    -R ${params.ref_fa} \
    -V ${vcf} \
    -O ${sampleID}_${output_suffix} \
    --cluster-window-size 10 \
    --filter-name "LowCoverage" --filter-expression "DP < 25" \
    --filter-name "LowQual" --filter-expression "QUAL < 20.0"
    """
}

// https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
