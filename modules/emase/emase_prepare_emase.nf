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

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-0f38b1b'

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
    emase prepare -G ${params.genome_file_list} -g ${params.gtf_file_list} -s ${params.haplotype_list} -o ./ -m --no-bowtie-index
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

 Usage: emase prepare [OPTIONS]

 prepare EMASE

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --genome-file      -G      FILE     Genome files, can seperate files by "," or have multiple -G [default: None] [required]                                                                                                                                      │
│    --haplotype-char   -s      TEXT     haplotype, either one per -h option, i.e. -h A -h B -h C, or a shortcut -h A,B,C [default: None]                                                                                                                            │
│    --gtf-file         -g      FILE     Gene Annotation File files, can seperate files by "," or have multiple -G [default: None]                                                                                                                                   │
│    --out_dir          -o      TEXT     Output folder to store results (default: the current working directory) [default: None]                                                                                                                                     │
│    --save-g2tmap      -m               saves gene id to transcript id list in a tab-delimited text file                                                                                                                                                            │
│    --save-dbs         -d               save dbs                                                                                                                                                                                                                    │
│    --no-bowtie-index  -m               skips building bowtie index                                                                                                                                                                                                 │
│    --verbose          -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                 │
│    --help                              Show this message and exit.                                                                                                                                                                                                 │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/