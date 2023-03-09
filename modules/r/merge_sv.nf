process MERGE_SV {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    tuple val(sampleID), val(meta), val(normal_name), val(tumor_name), file(manta_vcf), file(manta_vcf_tbi), val(manta), file(gridss_bgz), val(no_idx), val(gridss)
    val(chrom_list)

  output:
    tuple val(sampleID), file("${sampleID}.manta_gridss_sv.bed"), val(meta), emit: merged
    tuple val(sampleID), file("${sampleID}.manta_gridss_sv_supplemental.bed"), val(meta), emit: merged_suppl
    

  script:
    listOfChroms = chrom_list.collect { "$it" }.join(',')

    """
    Rscript ${projectDir}/bin/sv/merge-caller-vcfs.r \
        --vcf=${manta_vcf},${gridss_bgz} \
        --caller=manta,gridss \
        --tumor=${tumor_name} \
        --normal=${normal_name} \
        --build=GRCh38 \
        --slop=300 \
        --allowed_chr=${listOfChroms} \
        --min_sv_length=500 \
        --out_file=${sampleID}.manta_gridss_sv.bed \
        --out_file_supplemental=${sampleID}.manta_gridss_sv_supplemental.bed
    """
}