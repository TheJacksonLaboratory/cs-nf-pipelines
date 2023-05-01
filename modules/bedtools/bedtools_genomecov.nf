process BEDTOOLS_GENOMECOV {
  tag "$sampleID"

  cpus 2
  memory 4.GB 
  time '04:00:00'

  publishDir {
      def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control' : 'ip') : '' 
      "${params.pubdir}/${ params.organize_by=='sample' ? type+'/'+sampleID+'/bigwig' : 'bedtools'}"
  }, pattern: "*.txt", mode: 'copy'

 
  container 'quay.io/jaxcompsci/bedtools-sv_refs:2.30.0--hc088bd4_0'
 

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(flagstat)

  output:
  tuple val(sampleID), file("*.bedGraph"), emit: bedgraph
  tuple val(sampleID), file("*.txt"), emit: scale_factor


  script:
  log.info "----- Running bedtools genome coverage  on ${sampleID} -----"
  pe_fragment = params.read_type == 'SE' ? '' : '-pc'
  extend = (params.read_type == 'SE' && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
  """
  SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 1000000/\$1}')
  echo \$SCALE_FACTOR > ${sampleID}.scale_factor.txt
    
  bedtools genomecov -ibam ${bam[0]} -bg -scale \$SCALE_FACTOR $pe_fragment $extend | sort -T '.' -k1,1 -k2,2n > ${sampleID}.bedGraph
  """

}
