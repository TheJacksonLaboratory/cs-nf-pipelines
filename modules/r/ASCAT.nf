process ASCAT {
    tag "$sampleID"

    cpus 1
    memory 24.GB
    time '01:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}


    container 'quay.io/biocontainers/ascat:3.1.1--r43hdfd78af_1'

    input:
         val(sampleID),path(LRR),path(BAF),val(gender),val(platform),path(GC_file),path(RT_file)

    output:
        tuple val(sampleID),
             path("${sampleID}_sample.QC.txt"),
             path("${sampleID}_ASCAT_objects.Rdata"),
             path("${sampleID}.segments_raw.txt"),
             path("${sampleID}.segments.txt"),
             path("${sampleID}.aberrantcellfraction.txt"),
             path("${sampleID}.ploidy.txt"),
             path("ASCAT.failedarrays.txt", optional: true),
             path("ASCAT.nonaberrantarrays.txt", optional: true), emit: ascat

    script:
        """
        Rscript ${projectDir}/bin/cnv_array/ASCAT_run.R \
            ${sampleID} \
            ${LRR} \
            ${BAF} \
            ${meta.gender} \
            ${params.snp_platform} \
            ${params.GC_file} \
            ${params.RT_file}
        """
}