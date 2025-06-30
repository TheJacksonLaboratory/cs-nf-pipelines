process RSEM_PREPAREREFERENCE {
    cpus 12
    memory 64.GB
    time 12.h

    tag "$fasta $aligner $read_length"
    
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/rsem_bowtie2_star:0.2.0"

    publishDir "${outdir}", mode:'copy'

    input:
        tuple path(fasta), path(gtf), val(aligner), val(read_length)

    output:
        path("*"), emit: all_files
        path("*.transcripts.fa"), emit: transcripts

    script:
        mem_bytes = (task.memory.toMega() * 1024 * 1024) - 1000000 // Set star-ratelimit to module mem, and reserve 1GB for RSEM overhead.
        aligner_flags = (aligner == "bowtie2") ? "--bowtie2" : "--star --star-sjdboverhang ${read_length - 1} --star-ramlimit ${mem_bytes}"
        printf_statement = (aligner == "bowtie2") ? "'Prepared RSEM index for ${aligner} with ${fasta} and ${gtf}\n'" : "'Prepared RSEM index for ${aligner} with ${fasta} and ${gtf} and read length ${read_length}\n'"
        readme_filename = (aligner == "bowtie2") ? "README_${aligner}.txt" : "README_${aligner}_${read_length}.txt"
        outdir = (aligner == "bowtie2") ? "${params.pubdir}/${aligner}" : "${params.pubdir}/${aligner.toUpperCase()}/${fasta.baseName}_${read_length}"
        """
        rsem-prepare-reference \
            -p $task.cpus \
            --gtf ${gtf} \
            ${aligner_flags} \
            ${fasta} \
            ${fasta.baseName}

        printf ${printf_statement} > ${readme_filename}
        """
}

// Module adapted from: https://github.com/KU-GDSC/workflows
