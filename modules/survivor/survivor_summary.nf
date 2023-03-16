process SURVIVOR_SUMMARY {

    tag "$sampleID"

    cpus 1
    memory 2.GB
    time "00:30:00"

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    input:
        tuple val(sampleID), file(vcf)
    output:
        tuple val(sampleID), file("${sampleID}.survivor_summary.csv"), emit: csv
    script:
        """
        /usr/bin/env python ${projectDir}/bin/sv_to_table.py -v ${vcf} -o ${sampleID}.survivor_summary.csv
        """
}