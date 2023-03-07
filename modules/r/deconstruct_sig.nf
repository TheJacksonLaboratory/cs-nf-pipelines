process DECONSTRUCT_SIG {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '04:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bedpe'}/sigs", pattern: "*.txt", mode: 'copy'
  
  container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

  input:
    // COMPRESS_INDEX_MERGED_VCF.out.compressed_vcf_tbi
    tuple val(sampleID), file(vcf), file(idx), val(meta), val(normal_name), val(tumor_name)

  output:
    tuple val(sampleID), file("${sampleID}.cosmic.v3.2.deconstructSigs.signatures.highconfidence.txt"), emit: sigs
    tuple val(sampleID), file("${sampleID}.cosmic.v3.2.deconstructSigs.signatures.highconfidence.counts.txt"), emit: counts
    tuple val(sampleID), file("${sampleID}.cosmic.v3.2.deconstructSigs.signatures.highconfidence.input.txt"), emit: sigInput
    tuple val(sampleID), file("${sampleID}.cosmic.v3.2.deconstructSigs.signatures.highconfidence.reconstructed.txt"), emit: reconstructed
    tuple val(sampleID), file("${sampleID}.cosmic.v3.2.deconstructSigs.signatures.highconfidence.diff.txt"), emit: diff
  script:
    """
    Rscript ${projectDir}/bin/sv/run_deconstructSigs.R \
        --highconf "TRUE" \
        --file ${vcf} \
        --ref "GRCh38" \
        --cosmic ${cosmicSigs} \
        --output ${sampleID}.cosmic.v3.2.deconstructSigs.signatures.highconfidence \
        --samplename ${sampleID}
    """
}