process GBRS_RECONSTRUCT  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-16c7011'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.npz", mode: 'copy', enabled: params.keep_intermediate
    
    input:
    tuple val(sampleID), path(tpm), val(sex), val(generation)

    output:
    tuple val(sampleID), file("*genoprobs.npz"), emit: genoprobs_npz
    tuple val(sampleID), file("*genotypes.npz"), emit: genotypes_npz
    tuple val(sampleID), file("*genotypes.tsv"), emit: genotypes_tsv

    script:

    """
    cp ${params.base_ref_index_fai} ref.fa.fai

    gbrs reconstruct \
        -e ${tpm} \
        -t ${params.trans_prob_dir}/tranprob.DO.G${generation}.${sex}.npz \
        -x ${params.emission_prob_avecs} \
        -g ${params.gene_position_file} \
        -o ${sampleID} \
        -c ${params.gbrs_expression_threshold} \
        -s ${params.gbrs_sigma} 
    """

    stub:
    """
    touch ${sampleID}.genoprobs.npz
    touch ${sampleID}.genotypes.npz
    touch ${sampleID}.genotypes.tsv
    """
}

/*
 Usage: gbrs reconstruct [OPTIONS]

 reconstruct the genome based upon gene-level TPM quantities

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --expr-file       -e      FILE     file containing gene-level TPM quantities [default: None] [required]                                                                                                                                                         │
│ *  --tprob-file      -t      FILE     transition probabilities file [default: None] [required]                                                                                                                                                                     │
│    --avec-file       -x      FILE     alignment specificity file [default: None]                                                                                                                                                                                   │
│    --gpos-file       -g      FILE     meta information for genes (chrom, id, location) [default: None]                                                                                                                                                             │
│    --expr-threshold  -c      FLOAT    [default: 1.5]                                                                                                                                                                                                               │
│    --sigma           -s      FLOAT    [default: 0.12]                                                                                                                                                                                                              │
│    --outbase         -o      TEXT     basename of all the generated output files [default: None]                                                                                                                                                                   │
│    --verbose         -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                  │
│    --help                             Show this message and exit.                                                                                                                                                                                                  │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

 sigma defaults to 0.12.
 expr_threshold defaults to 1.5. 

*/
