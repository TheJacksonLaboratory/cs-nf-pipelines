process SAMTOOLS_FILTER_UNIQUE {
    tag "$sampleID"

    cpus  1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bicseq2:v2'

    input:
    tuple val(sampleID), val(meta), path(bam), path(bai), val(read_ID)
    val(chroms)

    output:
    tuple val(sampleID), path("seq_out/*.seq"), val(meta), val(read_ID), emit: uniq_seq

    script:
    chrom_list = chroms.collect { "$it" }.join(' ')
    """
    /samtools-0.1.7a_getUnique-0.1.3/samtools view -U "BWA,${read_ID}_,N,N" ${bam}
    
    mkdir seq_out

    for chrom in ${chrom_list}; do mv *_\$chrom.seq seq_out/; done
    """
}
// Modified samtools view: 
//          -U STR   If specified, get the uique reads. STR should be <Aligner,OutputPrefix,ChromNameReport?,StrandReport?> or <Aligner,OutputPrefix,ChromNameReport?,StrandReport?minLen,maxLen>
//                    e.g. <BWA,output,N,N> means that the aligner is BWA, the prefix of the output is output, and that the chromosome names and strands will not be reported
//                    StrandReport should be S (separate positive and negative reads), Y, or N
//                    minLen and maxLen specifies the length range of the reported reads

// NOTE: The modified samtools view with -U does not work as normal samtools. It does not pass the region to the filter step. 

// NOTE: The modified samtools view reads the header and uses the chrom names in the header, not matter what is present in the mapped file. 

// NOTE: Therefore, to avoid `chrUn_*_decoy.seq` and `_HLA-*.seq` non-primary chroms, the primary list is moved to an output directory. 
