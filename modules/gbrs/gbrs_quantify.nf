process GBRS_QUANTIFY {
    tag "$sampleID"

    cpus 1
    memory {5.GB * task.attempt}
    time {5.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.multiway.isoforms.tpm", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.multiway.isoforms.expected_read_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.multiway.isoforms.alignment_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.multiway.genes.tpm", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.multiway.genes.expected_read_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.multiway.genes.alignment_counts", mode: 'copy'

    input:
    tuple val(sampleID), path(h5)

    output:
    tuple val(sampleID), file("*.multiway.isoforms.tpm"), emit: isoforms_tpm
    tuple val(sampleID), file("*.multiway.isoforms.expected_read_counts"), emit: isoforms_expected_count
    tuple val(sampleID), file("*.multiway.isoforms.alignment_counts"), emit: isoforms_alignment_count
    tuple val(sampleID), file("*.multiway.genes.tpm"), emit: genes_tpm
    tuple val(sampleID), file("*.multiway.genes.expected_read_counts"), emit: genes_expected_cout
    tuple val(sampleID), file("*.multiway.genes.alignment_counts"), emit: genes_alignment_count

    script:
    """
    gbrs quantify \
        -i ${h5} \
        -g ${params.gene2transcript_csv} \
        -L ${params.full_transcript_info} \
        -M ${params.emase_model} \
        --report-alignment-counts \
        -o ${sampleID}
    """

    stub:
    """
    touch ${sampleID}.multiway.isoforms.tpm
    touch ${sampleID}.multiway.isoforms.expected_read_counts
    touch ${sampleID}.multiway.isoforms.alignment_counts
    touch ${sampleID}.multiway.genes.tpm
    touch ${sampleID}.multiway.genes.expected_read_counts
    touch ${sampleID}.multiway.genes.alignment_counts
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
