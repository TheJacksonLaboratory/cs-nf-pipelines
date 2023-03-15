process MANTA_CALL {
    tag "$sampleID"
    
    cpus = 8
    memory = 40.GB
    time = "10:00:00"

    container 'quay.io/biocontainers/manta:1.6.0--py27_0'

    publishDir "${params.outdir}/unmerged_calls", pattern: "${sampleID}_${caller}Sort.vcf", mode: 'copy'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
        tuple file(fasta), file(fai)
    
    output:
        tuple val(sampleID), file("${sampleID}_mantaSort.vcf"), emit: manta_sv

    script:
        """
        python /usr/local/bin/configManta.py \
            --runDir mantaSVOut \
            --bam ${bam} \
            --referenceFasta ${fasta}
        python ./mantaSVOut/runWorkflow.py -m local -j ${task.cpus}
        mv mantaSVOut/results/variants/candidateSV.vcf.gz ./${sampleID}_mantaCandidate.vcf.gz
        gunzip ${sampleID}_mantaCandidate.vcf.gz
        bash ${projectDir}/bin/vcfSort.sh ${sampleID}_mantaCandidate.vcf ${sampleID}_mantaSort.vcf
        """
}