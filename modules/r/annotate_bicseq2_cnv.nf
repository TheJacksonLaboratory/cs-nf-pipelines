process ANNOTATE_BICSEQ2_CNV {
  tag "$sampleID"

  cpus 1
  memory 36.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    tuple val(sampleID), file(bicseq2_calls), val(no_idx), val(meta), val(normal_name), val(tumor_name), val(bicseq2)
    val(chrom_list)

  output:
    tuple val(sampleID), file("${sampleID}.bicseq2_annotated.bed"), val(meta), val(normal_name), val(tumor_name), emit: bicseq_annot
    tuple val(sampleID), file("${sampleID}.bicseq2_annotated_supplemental.bed"), val(meta), val(normal_name), val(tumor_name), emit: bicseq_annot_suppl

  script:
    listOfChroms = chrom_list.collect { "$it" }.join(',')

    """
    Rscript ${projectDir}/bin/sv/annotate-cnv.r \
        --cnv = ${bicseq2_calls} \
        --caller = bicseq2 \
        --tumor = ${tumor_name} \
        --normal = ${normal_name} \
        --cytoband = ${params.cytoband} \
        --db_names=DGV,1000G,COSMIC \
        --db_files = ${params.dgv},${params.thousandG},${params.cosmicUniqueBed} \
        --cancer_census= = ${params.cancerCensusBed} \
        --ensembl = ${params.ensemblUniqueBed} \
        --allowed_chr = ${listOfChroms} \
        --overlap_fraction=0.8 \
        --out_file_main= ${sampleID}.bicseq2_annotated.bed \
        --out_file_supplemental = ${sampleID}.bicseq2_annotated_supplemental.bed

    """
}