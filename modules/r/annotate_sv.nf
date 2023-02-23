process ANNOTATE_SV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // MERGE_SV.out.merged
    tuple val(sampleID), file(merged_sv_bed), val(meta), val(normal_name), val(tumor_name)

  output:
    tuple val(sampleID), file("${sampleID}.manta_gridss_sv_annotated.bed"), val(meta), val(normal_name), val(tumor_name), emit: annot_sv_bed

  script:

    """
    Rscript /annotate-bedpe-with-databases.r \
        --db_names=gap,DGV,1000G,PON,COSMIC \
        --db_files=${params.gap},${params.dgvBedpe},${params.thousandGVcf},${params.svPon},${params.cosmicBedPe} \
        --slop=500 \
        --db_ignore_strand=COSMIC \
        --bedpe=${merged_sv_bed} \
        --out_file=${sampleID}.manta_gridss_sv_annotated.bed

    """
}