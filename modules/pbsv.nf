process DISCOVERTANDEM {

    tag "$sampleID"

    cpus 8
    memory { 40.GB * task.attempt }
    time { 4.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--h9ee0642_0'

    input:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai)
    output:
        tuple val(sampleID), file("${bam.baseName}.svsig.gz"), emit: pbsv_svsig
    script:
        """
        /usr/bin/env bash ${projectDir}/bin/pbsv_tandem.sh ${pbmm2_bam} ${pbmm2_bam.baseName}.svsig.gz
        """
}

process DISCOVERNOTANDEM {

    tag "$sampleID"

    cpus 8
    memory { 80.GB * task.attempt }
    time { 12.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--h9ee0642_0'

    input:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai)
    output:
        tuple val(sampleID), file("${bam.baseName}.svsig.gz"), emit: pbsv_svsig
    
    script:
        """
        pbsv discover ${pbmm2_bam} ${sampleID}.svsig.gz
        """
}

process CALLCCS {

    tag "$sampleID"

    cpus 8
    memory { 80.GB * task.attempt }
    time { 12.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--h9ee0642_0'

    publishDir "${params.pubdir}/unmerged_calls", pattern: "${sampleID}.pbsv_calls.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(pbsv_svsig)
        path(fasta)
    output:
        tuple val(sampleID), file("${sampleID}.pbsv_calls.vcf"), emit: pbsv_vcf
    script:
        """
        pbsv call --ccs ${fasta} ${pbsv_svsig} ${sampleID}.pbsv_calls.vcf
        """
}

process CALLCLR {

    tag "$sampleID"
    
    cpus 8
    memory { 80.GB * task.attempt }
    time { 12.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--h9ee0642_0'

    publishDir "${params.pubdir}/unmerged_calls", pattern: ${sampleID}.pbsv_calls.vcf, mode: "copy"

    input:
        tuple val(sampleID), file(pbsv_svsig)
        path(fasta)
    output:
        tuple val(sampleID), file("${sampleID}.pbsv_calls.vcf"), emit: pbsv_vcf
    script:
        """
        pbsv call ${fasta} ${pbsv_svsig} ${sampleID}.pbsv_calls.vcf
        """
}