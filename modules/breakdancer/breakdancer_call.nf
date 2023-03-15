process BREAKDANCER_CALL {
    tag "$sampleID"
    
    cpus = 1
    memory = 40.GB
    time = "10:00:00"

    container 'quay.io/biocontainers/breakdancer:1.4.5--h25a10a7_7'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
    
    output:
        tuple val(sampleID), file("${sampleID}_BreakDancer-SV"), emit: breakdancer_sv

    script:
        """
        bam2cfg.pl ${bam} > ${sampleID}_config
        breakdancer-max -r 5 -s 50 -h ${sampleID}_config > ${sampleID}_BreakDancer-SV
        """
}