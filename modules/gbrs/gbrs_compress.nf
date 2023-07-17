process GBRS_COMPRESS {
    tag "$sampleID"

    cpus 1
    memory { suffix == 'merged' ? 6.GB * task.attempt : 40.GB * task.attempt}
    time 10.hour
    errorStrategy 'finish' 

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-547132f'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.h5", mode: 'copy', enabled: "${ suffix == 'merged' || params.read_type == 'SE' ? true : false }"

    input:
    tuple val(sampleID), path(bam)
    val(suffix)

    output:
    tuple val(sampleID), file("*.compressed.emase.h5"), emit: compressed_emase_h5

    script:
    bam_list = bam.collect { "$it" }.join(' -i ')

    output_name = suffix == 'merged' ? "${sampleID}.merged.compressed.emase.h5" : "${bam[0].baseName}.compressed.emase.h5"

    """
    gbrs compress -i ${bam_list} -o ${output_name}
    """

    stub:
    """
    touch ${output_name}
    """
}

/*
 Usage: gbrs compress [OPTIONS]

 compress EMASE format alignment incidence matrix

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --emase-file  -i      FILE     EMASE file to compress, can seperate files by "," or have multiple -i [default: None] [required]                                                                                                                                 │
│ *  --output      -o      FILE     name of the compressed EMASE file [default: None] [required]                                                                                                                                                                     │
│    --comp-lib    -c      TEXT     compression library to use [default: zlib]                                                                                                                                                                                       │
│    --verbose     -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                      │
│    --help                         Show this message and exit.                                                                                                                                                                                                      │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/
