process GENERATE_GRID_FILE {

    cpus 2
    memory 40.GB
    time 48.hour

    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}", pattern: 'ref.genome_grid.GRCm39.tsv', mode:'copy'

    output:
    path('ref.genome_grid.GRCm39.tsv'), emit: h5_files

    script:
    """
    Rscript ${projectDir}/bin/gbrs/generate_grid_file.R
    """

    stub:
    """
    touch "ref.genome_grid.GRCm39.tsv"
    """

}
