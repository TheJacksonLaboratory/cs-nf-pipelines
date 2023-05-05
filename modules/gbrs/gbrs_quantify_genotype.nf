process GBRS_QUANTIFY_GENOTYPES {
    tag "$sampleID"

    cpus 1
    memory 5.GB
    time 5.hour
    errorStrategy 'finish' 
    maxRetries 1

    container 'quay.io/mikewlloyd/gbrs_test:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.diploid.isoforms.tpm", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.diploid.isoforms.expected_read_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.diploid.isoforms.alignment_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.diploid.genes.tpm", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.diploid.genes.expected_read_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.diploid.genes.alignment_counts", mode: 'copy'

    input:
    tuple val(sampleID), path(h5), path(genotype_tsv)

    output:
    tuple val(sampleID), file("*.diploid.isoforms.tpm"), emit: isoforms_tpm
    tuple val(sampleID), file("*.diploid.isoforms.expected_read_counts"), emit: isoforms_expected_count
    tuple val(sampleID), file("*.diploid.isoforms.alignment_counts"), emit: isoforms_alignment_count
    tuple val(sampleID), file("*.diploid.genes.tpm"), emit: genes_tpm
    tuple val(sampleID), file("*.diploid.genes.expected_read_counts"), emit: genes_expected_cout
    tuple val(sampleID), file("*.diploid.genes.alignment_counts"), emit: genes_alignment_count

    script:
    """
    gbrs quantify \
        -i ${h5} \
        -g ${params.gene2transcript_csv} \
        -G ${genotype_tsv} \
        -L ${params.full_transcript_info} \
        -M ${params.emase_model} \
        --report-alignment-counts \
        -o ${sampleID}
    """

    stub:
    """
    touch ${sampleID}.diploid.isoforms.tpm
    touch ${sampleID}.diploid.isoforms.expected_read_counts
    touch ${sampleID}.diploid.isoforms.alignment_counts
    touch ${sampleID}.diploid.genes.tpm
    touch ${sampleID}.diploid.genes.expected_read_counts
    touch ${sampleID}.diploid.genes.alignment_counts
    """
}

/*
NOTE: gbrs quantify is a wrapper around the `run-emase` code.

usage: gbrs quantify [-h] -i ALNFILE -g GRPFILE [-L LENFILE] [-G GTYPEFILE]
                     [-M MULTIREAD_MODEL] [-o OUTBASE] [-p PSEUDOCOUNT]
                     [-m MAX_ITERS] [-t TOLERANCE] [-a] [-w]

optional arguments:
  -h, --help            show this help message and exit
  -i ALNFILE, --alignment-file ALNFILE
  -g GRPFILE, --group-file GRPFILE
  -L LENFILE, --length-file LENFILE
  -G GTYPEFILE, --genotype GTYPEFILE
  -M MULTIREAD_MODEL, --multiread-model MULTIREAD_MODEL
  -o OUTBASE
  -p PSEUDOCOUNT, --pseudocount PSEUDOCOUNT
  -m MAX_ITERS, --max-iters MAX_ITERS
  -t TOLERANCE, --tolerance TOLERANCE
  -a, --report-alignment-counts
  -w, --report-posterior
*/
