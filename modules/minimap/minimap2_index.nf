process MINIMAP2_INDEX {

    tag "${fasta.baseName}"

    cpus 8
    memory 40.GB
    time "2:00:00"

    container 'quay.io-biocontainers-minimap2-2.24--h7132678_1'

    input:
        path(fasta)
    output:
        path "${fasta.baseName}.mmi", emit: minimap2_index
    script:
        """
        minimap2 -d ${fasta.baseName}.mmi ${fasta} 
        """
}