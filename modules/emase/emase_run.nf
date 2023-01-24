process EMASE_RUN {
    tag "$sampleID"

    cpus 1
    memory {5.GB * task.attempt}
    time {5.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.isoforms.tpm", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.isoforms.expected_read_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.isoforms.alignment_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.genes.tpm", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.genes.expected_read_counts", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.genes.alignment_counts", mode: 'copy'

    input:
    tuple val(sampleID), path(h5)

    output:
    tuple val(sampleID), file("*.isoforms.tpm"), emit: isoforms_tpm
    tuple val(sampleID), file("*.isoforms.expected_read_counts"), emit: isoforms_expected_count
    tuple val(sampleID), file("*.isoforms.alignment_counts"), emit: isoforms_alignment_count
    tuple val(sampleID), file("*.genes.tpm"), emit: genes_tpm
    tuple val(sampleID), file("*.genes.expected_read_counts"), emit: genes_expected_cout
    tuple val(sampleID), file("*.genes.alignment_counts"), emit: genes_alignment_count

    script:
    """
    run-emase \
        -i ${h5} \
        -g ${params.gene2transcript_list} \
        -L ${params.full_transcript_list} \
        -M ${params.emase_model} \
        -o ${sampleID}
    """

    stub:
    """
    touch ${sampleID}.isoforms.tpm
    touch ${sampleID}.isoforms.expected_read_counts
    touch ${sampleID}.isoforms.alignment_counts
    touch ${sampleID}.genes.tpm
    touch ${sampleID}.genes.expected_read_counts
    touch ${sampleID}.genes.alignment_counts
    """
}




/*
run-emase:
Usage:
    run-emase -i <h5_file> -g <grp_file> -L <len_file> -M <multiread_model> -o <outbase> \
              -p <pseudocount> -r <read_length> -m <max_iters> -t <tolerance>
Input:
    -i <h5_file>         : Alignments stored in a PyTables HDF5 format
    -g <grp_file>        : Gene-to-transcript map (ENSMUSGxxx followed by a list of ENSMUSTyyy's)
    -L <len_file>        : File that contains transcript lengths
    -M <multiread_model> : Multiread model ID
                           1: Gene->Allele->Isoform,
                           2: Gene->Isoform->Allele,
                           3: Gene->(Isoform*Allele),
                           4: Gene*Isoform*Allele  (default model)
    -o <outbase>         : EMASE outputs the result to <folder/basename> (default: './emase')
    -p <pseudocount>     : Pseudocount for allele specificity (default: 0.0)
    -r <read_length>     : Read length (default: 100)
    -m <max_iters>       : The number of maximum iterations for EM (default: 999)
    -t <tolerance>       : Tolerance for the termination of EM. (default: 0.0001)
Parameters:
    -h, --help : shows this help message
    -c         : reports the alignment counts (Consider using another script 'count-alignments' instead.)
    -w         : reports the posterior probability for each read
*/