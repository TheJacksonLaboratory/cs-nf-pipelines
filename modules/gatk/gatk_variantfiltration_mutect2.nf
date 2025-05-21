process GATK_VARIANTFILTRATION {
    tag "$sampleID"

    cpus = 1
    memory = 10.GB
    time = '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(vcf), file(idx)
    val(indel_snp)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf
    tuple val(sampleID), file("*.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    if (indel_snp == 'INDEL'){
        fs='200.0'
        output_suffix = 'INDEL_filtered.vcf'
    }
    if (indel_snp =='SNP'){
        fs ='60.0'
        output_suffix = 'SNP_filtered.vcf'
    }
    if (indel_snp == 'BOTH'){
        fs = '60.0'
        output_suffix = 'snp_indel_filtered.vcf'
    }

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" VariantFiltration \
    -R ${params.ref_fa} \
    -V ${vcf} \
    -O ${sampleID}_${output_suffix} \
    --cluster-window-size 10 \
    --filter-name "LowCoverage" --filter-expression "DP < 25" \
    --filter-name "StrandBias" --filter-expression "FS > ${fs}"
    """
}

// https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
