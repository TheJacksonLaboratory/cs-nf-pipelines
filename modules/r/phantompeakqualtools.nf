process PHANTOMPEAKQUALTOOLS {
    tag "$sampleID"
    
    cpus 8
    memory 50.GB
    time '04:00:00'

    container 'quay.io/biocontainers/phantompeakqualtools:1.2.2--0'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.out")  , emit: spp
    tuple val(sampleID), file("*.pdf")  , emit: pdf
    tuple val(sampleID), file("*.Rdata"), emit: rdata

    script:
    """
    RUN_SPP=`which run_spp.R`
    Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="$bam" -savp="${sampleID}.spp.pdf" -savd="${sampleID}.spp.Rdata" -out="${sampleID}.spp.out" -p=$task.cpus
    """
}
