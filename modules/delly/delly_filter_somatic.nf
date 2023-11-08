process DELLY_FILTER_SOMATIC {
    tag "$sampleID"
    
    cpus = 1
    memory 80.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/delly:1.1.6--h6b1aa3f_2'
    
    input:
    tuple val(sampleID), path(bcf), path(csi), val(meta), val(normal_name), val(tumor_name)
    
    output:
    tuple val(sampleID), path("*.bcf"), path("*.csi"), val(meta), val(normal_name), val(tumor_name), val('delly_sv'), emit: bcf_csi

    script:
    """
    echo -e "${normal_name}\tcontrol" > samples.tsv
    echo -e "${tumor_name}\ttumor" >> samples.tsv

    delly filter -f somatic -o ${sampleID}_delly_filtered_somaticSV.bcf -s samples.tsv ${bcf}
    """
}
