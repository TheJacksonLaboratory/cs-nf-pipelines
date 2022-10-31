process CONPAIR {
  tag "$sampleID"

  cpus 2
  memory 4.GB
  time '10:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'short_alignment_marking' }", pattern:"*.marked.bam", mode:'copy'
  
  input:
  tuple val(sampleID), file(tumor_bam)
  tuple val(sampleID), file(normal_bam)

  output:
  tuple val(sampleID), file("*.marked.bam"), emit: marked_bam

  script:
  log.info "----- Short Alignment Marking Running on: ${sampleID} -----"
  // estimate concordance and contamination estimator for tumor–normal pairs
  //Verifying concordance between two samples (tumor and normal)
  //Estimating contamination level in both the tumor and the normal:


  """
  /projects/compsci/jgeorge/harshpreet_CS/nextflow/pipeline/ngs-ops/nygc_pipeline/tools/Conpair/scripts/run_gatk_pileup_for_sample.py -B TUMOR_bam -O TUMOR_pileup --reference ${params.reference} --outfile ${sampleID}_pileup.txt
  /projects/compsci/jgeorge/harshpreet_CS/nextflow/pipeline/ngs-ops/nygc_pipeline/tools/Conpair/scripts/run_gatk_pileup_for_sample.py -B NORMAL_bam -O NORMAL_pileup --reference ${params.reference} --outfile ${sampleID}_pileup.txt

  /projects/compsci/jgeorge/harshpreet_CS/nextflow/pipeline/ngs-ops/nygc_pipeline/tools/Conpair/scripts/verify_concordance.py -T TUMOR_pileup -N NORMAL_pileup --outfile ${sampleID}_concordance.txt

  /projects/compsci/jgeorge/harshpreet_CS/nextflow/pipeline/ngs-ops/nygc_pipeline/tools/Conpair/scripts/estimate_tumor_normal_contamination.py -T TUMOR_pileup -N NORMAL_pileup --outfile ${sampleID}_contamination.txt

  """

}
