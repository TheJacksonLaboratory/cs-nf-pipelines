process EMASE_CREATE_HYBRID {
    
    // give a group of fastas, and a haplotype list 
    // 1. generate a hybrid genome
    // 2. generate transcript list
    // 3. generate bowtie index. 

    cpus 1
    memory 15.GB
    time 5.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-16c7011'

    publishDir "${params.pubdir}/emase", pattern: "*.fa", mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: "*.info", mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: "*.tsv", mode:'copy'
    publishDir "${params.pubdir}/emase/bowtie", pattern: "*.ebwt", mode:'copy'

    output:
    path file("*.fa"), emit: transcript_fasta
    path file("*.info"), emit: transcript_info
    path file("*.ebwt"), emit: bowtie_index

    script:
    """
    emase create-hybrid -F ${params.genome_file_list} -s ${params.haplotype_list} -o ./ --create-bowtie-index
    """

    stub:
    """
    touch emase.pooled.transcripts.fa
    touch emase.pooled.transcripts.info
    touch bowtie.transcripts.4.ebwt
    touch bowtie.transcripts.3.ebwt
    touch bowtie.transcripts.2.ebwt
    touch bowtie.transcripts.1.ebwt
    touch bowtie.transcripts.rev.2.ebwt
    touch bowtie.transcripts.rev.1.ebwt
    """
}


/*
 Usage: emase create-hybrid [OPTIONS]

 hybridize Fasta files

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --target-files        -F                             FILE     Fasta file to parse, can seperate files by "," or have multiple -i [default: None] [required]                                                                                                     │
│ *  --suffices            -s                             TEXT     haplotype, either one per -h option, i.e. -h A -h B -h C, or a shortcut -h A,B,C [default: None] [required]                                                                                       │
│    --output              -o                             FILE     [default: gbrs.hybridized.targets.fa]                                                                                                                                                             │
│    --build-bowtie-index      --no-build-bowtie-index             [default: no-build-bowtie-index]                                                                                                                                                                  │
│    --verbose             -v                             INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                       │
│    --help                                                        Show this message and exit.                                                                                                                                                                       │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/