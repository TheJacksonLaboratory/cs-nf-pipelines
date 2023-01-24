process GBRS_BAM2EMASE {
    tag "$sampleID"

    cpus 1
    memory {15.GB * task.attempt}
    time {5.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gbrs' }", pattern: "*.h5", mode: 'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.emase.h5"), emit: emase_h5

    script:
    """
    gbrs bam2emase -i ${bam} \
                -m ${params.transcripts_info} \
                -s ${params.gbrs_strain_list} \
                -o ${bam.baseName}.emase.h5
    """

    stub:
    """
    touch ${bam.baseName}.emase.h5
    """
}


/*
usage: gbrs bam2emase [-h] -i ALNFILE [-m LIDFILE] [-s HAPLOGYPES]
                      [-o OUTFILE] [--delim DELIM] [--index-dtype INDEX_DTYPE]
                      [--data-dtype DATA_DTYPE]

optional arguments:
  -h, --help            show this help message and exit
  -i ALNFILE, --bamfile ALNFILE
  -m LIDFILE, --locus-ids LIDFILE
  -s HAPLOGYPES, --haplotype-codes HAPLOGYPES
  -o OUTFILE
  --delim DELIM
  --index-dtype INDEX_DTYPE
  --data-dtype DATA_DTYPE
*/