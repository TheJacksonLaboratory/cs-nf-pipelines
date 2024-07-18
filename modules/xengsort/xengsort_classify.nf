process XENGSORT_CLASSIFY {
    
    tag "$sampleID"

    // resource utilization
    cpus 32
    memory 60.GB
    time 48.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // load xengsort container
    container 'quay.io/biocontainers/xengsort:2.0.5--pyhdfd78af_0'

    // output directory
    // publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/xengsort/xengsort_classify' : 'xengsort'}", pattern: "*.fq", mode: "copy"
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/stats': 'xengsort' }", pattern: "*.txt", mode:'copy'

    // inputs
    input:
    path(xengsort_index)
    tuple val(sampleID), path(trimmed)

    output:
    tuple val(sampleID), path("*fastq-graft_sorted.*.fq"), emit: xengsort_human_fastq
    tuple val(sampleID), path("*fastq-host_sorted.*.fq"), emit: xengsort_mouse_fastq
    tuple val(sampleID), path("*.txt"), emit: xengsort_log
    
    script:

    // specify single-end or paired-end data
    if (params.read_type == "SE")

        """
        
        xengsort classify \
        --index ${xengsort_index}/${xengsort_index} \
        --fastq ${trimmed[0]} \
        --prefix ${sampleID} \
        --mode count \
        --threads ${task.cpus}
        --out fastq \
        --chunksize 32.0 \
        --compression none &> ${sampleID}_xengsort_log.txt

        cat fastq-host.1.fq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sampleID}_fastq-host_sorted.1.fq
        cat fastq-graft.1.fq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sampleID}_fastq-graft_sorted.1.fq

        """

    else if (params.read_type == "PE")

        """

        xengsort classify \
        --index ${xengsort_index}/${params.xengsort_idx_name} \
        --fastq ${trimmed[0]} \
        --pairs ${trimmed[1]} \
        --prefix ${sampleID} \
        --mode count \
        --threads ${task.cpus} \
        --out fastq \
        --chunksize 32.0 \
        --compression none &> ${sampleID}_xengsort_log.txt

        cat fastq-host.1.fq | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" > ${sampleID}_fastq-host_sorted.1.fq
        cat fastq-host.2.fq | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" > ${sampleID}_fastq-host_sorted.2.fq

        cat fastq-graft.1.fq | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" > ${sampleID}_fastq-graft_sorted.1.fq
        cat fastq-graft.2.fq | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" > ${sampleID}_fastq-graft_sorted.2.fq

        """

    else error "${params.read_type} is invalid, specify either SE or PE"
}