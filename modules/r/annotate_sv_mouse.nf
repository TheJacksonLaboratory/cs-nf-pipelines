process ANNOTATE_SV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // MERGE_SV.out.merged
    tuple val(sampleID), file(merged_sv_bed), val(normal_name), val(tumor_name)
    val(suppl_switch)

  output:
    tuple val(sampleID), file("${sampleID}.manta_lumpy_delly_sv_annotated*.bed"), val(normal_name), val(tumor_name), emit: annot_sv_bedpe

  script:

    if (suppl_switch == "main")
    """
    Rscript ${projectDir}/bin/pta/annotate-bedpe-with-databases.r \
        --db_names=GAP,MGP,EXCLUDE_RANGE \
        --db_files=${params.gap},${params.mgp},${params.exclude_list} \
        --slop=500 \
        --bedpe=${merged_sv_bed} \
        --out_file=${sampleID}.manta_lumpy_delly_sv_annotated.bed

    """
    else if (suppl_switch == "supplemental")
    """
    Rscript ${projectDir}/bin/pta/annotate-bedpe-with-databases.r \
        --db_names=GAP,MGP,EXCLUDE_RANGE \
        --db_files=${params.gap},${params.mgp},${params.exclude_list} \
        --slop=500 \
        --bedpe=${merged_sv_bed} \
        --out_file=${sampleID}.manta_lumpy_delly_sv_annotated_supplemental.bed
    """
}