process ANNOTATE_SV_WITH_CNV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    tuple val(sampleID), val(meta), file(bicseq_annot), file(annot_sv_genes_bedpe)

  output:
    tuple val(sampleID), file("${sampleID}.manta_gridss_sv_annotated_genes_cnv.bed"), val(meta), emit: sv_genes_cnv_bedpe
 
  script:

    """
    Rscript ${projectDir}/bin/sv/annotate-bedpe-with-cnv.r \
        --cnv=${bicseq_annot} \
        --bedpe=${annot_sv_genes_bedpe} \
        --out_file=${sampleID}.manta_gridss_sv_annotated_genes_cnv.bed
    """
}