process FILTER_BEDPE {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bedpe'}", pattern: "*.bedpe", mode: 'copy'
  
  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // ANNOTATE_SV_WITH_CNV.out.sv_genes_cnv_bedpe
    tuple val(sampleID), file(sv_genes_cnv_bedpe), val(normal_name), val(tumor_name)
    val(suppl_switch)
  output:
    tuple val(sampleID), file("${sampleID}_sv_annotated_somatic_final.bedpe"), val(normal_name), val(tumor_name), optional: true
    tuple val(sampleID), file("${sampleID}_sv_annotated_somatic_supplemental.bedpe"), val(normal_name), val(tumor_name), optional: true
    tuple val(sampleID), file("${sampleID}_sv_annotated_somatic_high_confidence_final.bedpe"), val(normal_name), val(tumor_name), optional: true
    tuple val(sampleID), file("${sampleID}_sv_annotated_somatic_high_confidence_supplemental.bedpe"), val(normal_name), val(tumor_name), optional: true
 
  script:
    if(suppl_switch == "main")
    """
    Rscript ${projectDir}/bin/pta/filter-bedpe.r \
        --max_changepoint_distance=1000 \
        --filter_databases=INS,DEL,INV \
        --bedpe=${sv_genes_cnv_bedpe} \
        --genome=GRCm39 \
        --out_file_somatic=${sampleID}_sv_annotated_somatic_final.bedpe \
        --out_file_highconf=${sampleID}_sv_annotated_somatic_high_confidence_final.bedpe
    """

    else if (suppl_switch == "supplemental")
    """
    Rscript ${projectDir}/bin/pta/filter-bedpe.r \
        --max_changepoint_distance=1000 \
        --filter_databases=INS,DEL,INV \
        --bedpe=${sv_genes_cnv_bedpe} \
        --genome=GRCm39 \
        --out_file_somatic=${sampleID}_sv_annotated_somatic_supplemental.bedpe \
        --out_file_highconf=${sampleID}_sv_annotated_somatic_high_confidence_supplemental.bedpe
    """
}
