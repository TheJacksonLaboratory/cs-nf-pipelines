process GBRS_QUANTIFY_GENOTYPES {
    tag "$sampleID"

    cpus 1
    memory 5.GB
    time 5.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-547132f'

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
        -a \
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
