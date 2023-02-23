process ANNOTATE_SV_WITH_CNV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // ANNOTATE_BICSEQ2_CNV.out.bicseq_annot
    tuple val(sampleID), file(bicseq_annot), val(meta), val(normal_name), val(tumor_name)
    // ANNOTATE_GENES_SV.out.annot_sv_genes_bedpe
    tuple val(sampleID), file(annot_sv_genes_bedpe), val(meta), val(normal_name), val(tumor_name)

  output:
    tuple val(sampleID), file("${sampleID}.manta_gridss_sv_annotated_genes_cnv.bed"), val(meta), val(normal_name), val(tumor_name), emit: sv_genes_cnv_bedpe
 
  script:

    """
    Rscript ${projectDir}/bin/sv/annotate-bedpe-with-cnv.r \
        --cnv=${bicseq_annot} \
        --bedpe=${annot_sv_genes_bedpe} \
        --out_file=${sampleID}.manta_gridss_sv_annotated_genes_cnv.bed
    """
}