process PBMM2_INDEX {

    tag "${fasta.baseName}"

    cpus 8
    memory 40.GB
    time "2:00:00"

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        path(fasta)
    output:
        path "${fasta.baseName}.mmi", emit: pbmm2_index
    script:
        """
        pbmm2 index ${fasta} ${fasta.baseName}.mmi
        """
}