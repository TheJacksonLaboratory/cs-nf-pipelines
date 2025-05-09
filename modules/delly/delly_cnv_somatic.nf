process DELLY_CNV_SOMATIC {
    tag "$sampleID"
    
    cpus = 1
    memory 80.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/delly:1.1.6--h6b1aa3f_2'
    
    publishDir "${params.pubdir}/${sampleID + '/cnv_plots'}", pattern: "*.cov.gz", mode: 'copy'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), path("${normal_name}.bcf"), path("${tumor_name}.bcf"), val(meta), val(normal_name), val(tumor_name), emit: bcfs
    tuple val(sampleID), path("${tumor_name}.cov.gz"), emit: tumor_cov


    script:
    """
    delly cnv -u -z ${params.cnv_min_size} -i ${params.cnv_window} -o ${tumor_name}.bcf -c ${tumor_name}.cov.gz -g ${params.ref_fa} -m ${params.delly_mappability} ${tumor_bam} 
    delly cnv -u -v ${tumor_name}.bcf -o ${normal_name}.bcf -g ${params.ref_fa} -m ${params.delly_mappability} ${normal_bam}
    """

    stub:
    """
    wget -O ${tumor_name}.bcf https://raw.githubusercontent.com/TheJacksonLaboratory/cs-nf-test/main/pta/mouse/fizzbang--n_fizz.bcf
    wget -O ${normal_name}.bcf https://raw.githubusercontent.com/TheJacksonLaboratory/cs-nf-test/main/pta/mouse/fizzbang--t_bang.bcf
    wget -O ${tumor_name}.cov.gz https://raw.githubusercontent.com/TheJacksonLaboratory/cs-nf-test/main/pta/mouse/fizzbang--t_bang.cov.gz
    """
}

// Delly coverage requirements preculde most testing with small sample files. 
// The stub run is included here to overcome the coverage limitations. 
// Additional testing of the module has been done, and will be written formally in the future.
// note that due to sample naming, the stub run must be done with fizzbang--n_fizz and fizzbang--t_bang