process XENGSORT_CLASSIFY {
    
    tag "$sampleID"

    // resource utilization
    cpus 32
    memory 60.GB
    time 48.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // load xengsort container
    container 'quay.io/jaxcompsci/xengsort_gnu_utils:v2.0.5'

    // output directory
    // publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/xengsort/xengsort_classify' : 'xengsort'}", pattern: "*.fq", mode: "copy"
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/stats': 'xengsort' }", pattern: "*.txt", mode:'copy'

    // inputs
    input:
    path(xengsort_index)
    tuple val(sampleID), path(trimmed)

    output:
    tuple val(sampleID), path("*graft_sorted.*.fq"), emit: xengsort_human_fastq
    tuple val(sampleID), path("*host_sorted.*.fq"), emit: xengsort_mouse_fastq
    tuple val(sampleID), path("*.txt"), emit: xengsort_log
    
    script:

    // specify single-end or paired-end data
    if (params.read_type == "SE")

        """
        
        xengsort classify \
        --index ${xengsort_index}/${params.xengsort_idx_name} \
        --fastq ${trimmed[0]} \
        --prefix ${sampleID} \
        --mode count \
        --threads ${task.cpus} \
        --chunksize 32.0 \
        --compression none &> ${sampleID}_xengsort_log.txt

        cat ${sampleID}-host.fq | paste - - - - | sort -k1,1 -T ./ -t " " | tr "\\t" "\\n" > ${sampleID}-host_sorted.1.fq
        cat ${sampleID}-graft.fq | paste - - - - | sort -k1,1 -T ./ -t " " | tr "\\t" "\\n" > ${sampleID}-graft_sorted.1.fq

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
        --chunksize 32.0 \
        --compression none &> ${sampleID}_xengsort_log.txt

        cat ${sampleID}-host.1.fq | paste - - - - | sort -k1,1 -T ./ -t " " | tr "\\t" "\\n" > ${sampleID}-host_sorted.1.fq
        cat ${sampleID}-host.2.fq | paste - - - - | sort -k1,1 -T ./ -t " " | tr "\\t" "\\n" > ${sampleID}-host_sorted.2.fq

        cat ${sampleID}-graft.1.fq | paste - - - - | sort -k1,1 -T ./ -t " " | tr "\\t" "\\n" > ${sampleID}-graft_sorted.1.fq
        cat ${sampleID}-graft.2.fq | paste - - - - | sort -k1,1 -T ./ -t " " | tr "\\t" "\\n" > ${sampleID}-graft_sorted.2.fq

        """

    else error "${params.read_type} is invalid, specify either SE or PE"
}