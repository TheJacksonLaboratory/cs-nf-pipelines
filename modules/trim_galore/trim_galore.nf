process TRIM_GALORE {
    tag "$sampleID"

    cpus 8
    memory 16.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${type + sampleID + '/trimmed_fastq'}"
    }, pattern: "*.fq.gz", mode: 'copy', enabled: params.keep_intermediate

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? 'fastqc/' : ''
        "${params.pubdir}/${type + sampleID + '/stats'}"
    }, pattern: "*_fastqc.{zip,html}", mode: 'copy' 

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? 'fastqc/' : ''
        "${params.pubdir}/${type + sampleID + '/trimmed_fastq'}"
    }, pattern: "*trimming_report.txt", mode: 'copy'


    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*_fastqc.{zip,html}"), emit: trimmed_fastqc
    tuple val(sampleID), path("*.fq.gz"), emit: trimmed_fastq
    tuple val(sampleID), path("*trimming_report.txt"), emit: trim_stats

    script:

    paired_end = params.read_type == 'PE' ?  '--paired' : ''
    rrbs_flag = params.workflow == "rrbs" ? '--rrbs' : ''
    /*
    RRBS Mode
    In this mode, Trim Galore identifies sequences that were adapter-trimmed and removes another 2 bp from the 3' end of Read 1, 
    and for paired-end libraries also the first 2 bp of Read 2 (which is equally affected by the fill-in procedure). 
    This is to avoid that the filled-in cytosine position close to the second MspI site in a sequence is used for methylation calls. 
    Sequences which were merely trimmed because of poor quality will not be shortened any further.
    */
    directionality = params.non_directional ? '--non_directional': ''
    /*
    Non-directional mode
    The non-directional option, will screen adapter-trimmed sequences for the presence of either CAA or CGA at the start of sequences
    and clip off the first 2 bases if found. If CAA or CGA are found at the start, no bases will be trimmed off from the 3’ end 
    even if the sequence had some contaminating adapter sequence removed   in this case the sequence read likely originated from either the CTOT or CTOB strand; 
    refer to the RRBS guide for the meaning of CTOT and CTOB strands). 
    */

    if (params.workflow == "chipseq" && params.read_type == 'SE')
    """
        [ ! -f  ${sampleID}.fastq.gz ] && ln -s ${fq_reads} ${sampleID}.fastq.gz

        trim_galore --cores ${task.cpus} ${paired_end} ${rrbs_flag} ${directionality} --gzip --length ${params.trimLength} -q ${params.qualThreshold}  --stringency ${params.adapOverlap}  -a ${params.adaptorSeq}  --fastqc ${sampleID}.fastq.gz
    """
    else if (params.workflow == "chipseq" && params.read_type == 'PE')
    """
        [ ! -f  ${sampleID}_1.fastq.gz ] && ln -s ${fq_reads[0]} ${sampleID}_1.fastq.gz
        [ ! -f  ${sampleID}_2.fastq.gz ] && ln -s ${fq_reads[1]} ${sampleID}_2.fastq.gz

        trim_galore --cores ${task.cpus} ${paired_end} ${rrbs_flag} ${directionality} --gzip --length ${params.trimLength} -q ${params.qualThreshold}  --stringency ${params.adapOverlap}  -a ${params.adaptorSeq}  --fastqc ${sampleID}_1.fastq.gz ${sampleID}_2.fastq.gz
    """
    else
    """ 
        trim_galore --basename ${sampleID} --cores ${task.cpus} ${paired_end} ${rrbs_flag} ${directionality} --gzip --length ${params.trimLength} -q ${params.qualThreshold}  --stringency ${params.adapOverlap}  -a ${params.adaptorSeq}  --fastqc ${fq_reads}
    """
}
