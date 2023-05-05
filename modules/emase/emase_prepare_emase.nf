process EMASE_PREPARE_EMASE {

    // give a fasta or group of fastas, gtf or group of gtfs, and a haplotype list 
    // 1. generate a hybrid genome
    // 2. generate transcript list
    // 3. generate gene to transcript map
    // 4. generate bowtie index. 

    // NOTE: Transcript lists are 'pooled' but can be incomplete for certain haplotypes. 
    //       Missing transcripts in haplotypes will cause errors in `run-emase`. 
    //       Helper script `clean_transcript_info.py` can be used to add missing transcripts. 

    cpus 1
    memory 15.GB
    time 24.hour

    container 'quay.io/mikewlloyd/gbrs_test:latest'

    publishDir "${params.pubdir}/emase", pattern: '*.fa', mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: '*.info', mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/emase", pattern: '*.tsv', mode:'copy'
    // publishDir "${params.pubdir}/bowtie", pattern: "*.ebwt", mode:'copy' // TURN ON IF BOWTIE INDEX BUILT HERE.

    output:
    path("*.fa"), emit: pooled_transcript_fasta
    path("*.info"), emit: pooled_transcript_info
    path("*.tsv"), emit: pooled_gene_to_transcripts
    // path("*.ebwt"), emit: pooled_bowtie_index, optional: true // TURN ON IF BOWTIE INDEX BUILT HERE.

    script:
    """
    prepare-emase -G ${params.genome_file_list} -g ${params.gtf_file_list} -s ${params.haplotype_list} -o ./ -m --no-bowtie-index
    """

    stub:
    """
    touch emase.pooled.transcripts.fa
    touch emase.pooled.transcripts.info
    touch emase.gene2transcripts.tsv
    touch bowtie.transcripts.4.ebwt
    touch bowtie.transcripts.3.ebwt
    touch bowtie.transcripts.2.ebwt
    touch bowtie.transcripts.1.ebwt
    touch bowtie.transcripts.rev.2.ebwt
    touch bowtie.transcripts.rev.1.ebwt
    """

}

/*
prepare-emase:
Usage:
    prepare-emase -G <genome_files> [ -g <gtf_files> -s <hap_list> -o <out_dir> -m -x ]

Input:
    -G <genome_files> : List of Genome files (comma delimited)
    -g <gtf_files>    : List of gene annotation files (comma delimited, in the order of genomes)
    -s <hap_list>     : Names of haplotypes to be used instead (comma delimited, in the order of genomes)
    -o <out_dir>      : Output folder to store results (default: the current working directory)

Parameters:
    -h, --help            : shows this help message
    -m, --save-g2tmap     : saves gene id to transcript id list in a tab-delimited text file
    -x, --no-bowtie-index : skips building bowtie index

Note:
    Does not work if the input gtf file is older than Ensembl Release 75.
*/