process BISMARK_METHYLATION_EXTRACTION {
    tag "$sampleID"

    cpus 8
    memory 60.GB
    time 30.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bismark:0.23.1--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*txt", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/bismark_results'}", pattern: "*.{png,gz}", mode:'copy'  

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*splitting_report.txt"), emit: extractor_reports
    tuple val(sampleID), file("*.M-bias.txt"), emit: extractor_mbias
    tuple val(sampleID), file("*.{png,gz}"), emit: extractor_png_gz

    script:
    
    comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
    cytosine_report = params.cytosine_report ? "--cytosine_report --genome_folder ${params.ref_fa_index}" : ''
    
    if(params.read_type == 'SE') {
        """
        bismark_methylation_extractor  -multicore ${task.cpus} ${comprehensive} ${cytosine_report} --bedGraph --counts --gzip -s --report ${bam}
        """
    } else {
        """
        bismark_methylation_extractor -multicore ${task.cpus} ${comprehensive} ${cytosine_report} --ignore_r2 2 --ignore_3prime_r2 2 --bedGraph --counts --gzip -p --no_overlap --report ${bam}
        """
    }
}
