process FILTER_BEDPE {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bedpe'}", patterh: "*.bedpe", mode: 'copy'
  
  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // ANNOTATE_SV_WITH_CNV.out.sv_genes_cnv_bedpe
    tuple val(sampleID), file(sv_genes_cnv_bedpe), val(meta)
    val(suppl_switch)
  output:
    tuple val(sampleID), file("${sampleID}.sv.annotated.v7.somatic.final.bedpe"), val(meta), optional: true
    tuple val(sampleID), file("${sampleID}.sv.annotated.v7.somatic.supplemental.bedpe"), val(meta), optional: true
    tuple val(sampleID), file("${sampleID}.sv.annotated.v7.somatic.high_confidence.final.bedpe"), val(meta), optional: true
    tuple val(sampleID), file("${sampleID}.sv.annotated.v7.somatic.high_confidence.supplemental.bedpe"), val(meta), optional: true
 
  script:
    if(suppl_switch == "main")
    """
    Rscript ${projectDir}/bin/sv/filter-bedpe.r \
        --max_changepoint_distance=1000 \
        --filter_databases=DGV,1000G,PON \
        --bedpe=${sv_genes_cnv_bedpe} \
        --out_file_somatic=${sampleID}.sv.annotated.v7.somatic.final.bedpe \
        --out_file_highconf=${sampleID}.sv.annotated.v7.somatic.high_confidence.final.bedpe
    """

    else if (suppl_switch == "supplemental")
    """
    Rscript ${projectDir}/bin/sv/filter-bedpe.r \
        --max_changepoint_distance=1000 \
        --filter_databases=DGV,1000G,PON \
        --bedpe=${sv_genes_cnv_bedpe} \
        --out_file_somatic=${sampleID}.sv.annotated.v7.somatic.supplemental.bedpe \
        --out_file_highconf=${sampleID}.sv.annotated.v7.somatic.high_confidence.supplemental.bedpe
    """
}