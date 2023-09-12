process MERGE_SV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    tuple val(sampleID), val(normal_name), val(tumor_name), path(manta_vcf), path(manta_tbi), path(lumpy_vcf), path(lumpy_tbi), path(delly_vcf), path(delly_tbi)
    val(chrom_list)

  output:
    tuple val(sampleID), path("${sampleID}.manta_lumpy_delly_sv.bed"), val(normal_name), val(tumor_name), emit: merged
    tuple val(sampleID), path("${sampleID}.manta_lumpy_delly_sv_supplemental.bed"), val(normal_name), val(tumor_name), emit: merged_suppl
    

  script:
    listOfChroms = chrom_list.collect { "$it" }.join(',')

    """
    Rscript ${projectDir}/bin/pta/merge-caller-vcfs.r \
        --vcf=${manta_vcf},${lumpy_vcf},${delly_vcf} \
        --caller=manta,lumpy,delly \
        --tumor=${tumor_name} \
        --normal=${normal_name} \
        --build=GRCm39 \
        --slop=1000 \
        --allowed_chr=${listOfChroms} \
        --min_sv_length=200 \
        --out_file=${sampleID}.manta_lumpy_delly_sv.bed \
        --out_file_supplemental=${sampleID}.manta_lumpy_delly_sv_supplemental.bed
    """
}
