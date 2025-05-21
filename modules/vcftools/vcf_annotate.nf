process VCF_ANNOTATE {
    tag "$sampleID"

    cpus = 1
    memory = 10.GB
    time = '23:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    input:
    tuple val(sampleID), file(snp_vcf)
    val(indel_snp)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    container 'quay.io/biocontainers/perl-vcftools-vcf:0.1.16--pl5321hdfd78af_4'

    script:

    if (params.gen_org=='mouse'){
        delta="CHROM,POS,ID,REF,ALT"
    }
    else if (params.gen_org=='human'){
        delta="CHROM,POS,ID,REF,ALT"
    }

    if (indel_snp == 'INDEL'){
        output_suffix = 'indel_filtered_dbsnp.vcf'
    }
    if (indel_snp =='SNP'){
        output_suffix = 'snp_filtered_dbsnp.vcf'
    }
    if (indel_snp == 'BOTH'){
        output_suffix = 'snp_indel_filtered_dbsnp.vcf'
    }

    """
    cat ${snp_vcf} | vcf-annotate -a ${params.dbSNP} -c ${delta} > ${sampleID}_${output_suffix}
    """
}
