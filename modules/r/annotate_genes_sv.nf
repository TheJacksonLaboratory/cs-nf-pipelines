process ANNOTATE_GENES_SV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // ANNOTATE_SV.out.annot_sv_bedpe
    tuple val(sampleID), file(annot_sv_bedpe), val(meta), val(normal_name), val(tumor_name)

  output:
    tuple val(sampleID), file("${sampleID}.manta_gridss_sv_annotated_genes.bed"), val(meta), val(normal_name), val(tumor_name), emit: annot_sv_genes_bedpe

  script:

    """
    Rscript ${projectDir}/bin/sv/annotate-bedpe-with-genes.r \
        --ensembl=${params.ensemblUniqueBed} \
        --cancer_census=${params.cancerCensusBed} \
        --bedpe=${annot_sv_bedpe} \
        --out_file=${sampleID}.manta_gridss_sv_annotated_genes.bed

    """
}