process GBRS_QUANTIFY {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time 2.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-16c7011'

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
        -a \
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

 Usage: gbrs quantify [OPTIONS]

 quantify allele-specific expressions

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --alignment-file           -i      FILE     EMASE alignment incidence file (in hdf5 format) [default: None] [required]                                                                                                                                          │
│    --group-file               -g      FILE     tab delimited file of gene to transcript mapping [default: None]                                                                                                                                                    │
│    --length-file              -L      FILE     tab delimited file of locus(transcript) and length [default: None]                                                                                                                                                  │
│    --genotype                 -G      FILE     tab delimited file of locus(transcipt) and diplotype [default: None]                                                                                                                                                │
│    --outbase                  -o      TEXT     basename of all the generated output files [default: gbrs.quantified]                                                                                                                                               │
│    --multiread-model          -M      INTEGER  emase model (default: 4) [default: 4]                                                                                                                                                                               │
│    --pseudocount              -p      FLOAT    prior read count (default: 0.0) [default: 0.0]                                                                                                                                                                      │
│    --max-iters                -m      INTEGER  maximum iterations for EM iteration [default: 999]                                                                                                                                                                  │
│    --tolerance                -t      FLOAT    tolerance for EM termination (default: 0.0001 in TPM) [default: 0.0001]                                                                                                                                             │
│    --report-alignment-counts  -a               whether to report alignment counts                                                                                                                                                                                  │
│    --report-posterior         -w               whether to report posterior probabilities                                                                                                                                                                           │
│    --verbose                  -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                         │
│    --help                                      Show this message and exit.                                                                                                                                                                                         │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/
