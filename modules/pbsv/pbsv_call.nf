process PBSV_CALL {

    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--h9ee0642_0'

    publishDir "${params.pubdir}/unmerged_calls", pattern: "${sampleID}.pbsv_calls.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(pbsv_svsig)
        path(fasta)
    output:
        tuple val(sampleID), file("${sampleID}.pbsv_calls.vcf"), emit: pbsv_vcf
    script:
        if (params.pbmode == "CCS")
            """
            pbsv call --ccs ${fasta} ${pbsv_svsig} ${sampleID}.pbsv_calls.vcf
            """
        else if (params.pbmode == "CLR")
            """
            pbsv call ${fasta} ${pbsv_svsig} ${sampleID}.pbsv_calls.vcf
            """
}