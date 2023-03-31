process CALC_PBC_METRICS {
  tag "$sampleID"

  cpus 4
  memory 20.GB    
  time = '10:00:00' 

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'bedtools' }", pattern: "*.pbc.qc", mode: 'copy'
  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

  input:
  tuple val(sampleID), file(tmp_bams)

  output:
  tuple val(sampleID), file("*.pbc.qc")

  shell:
  log.info "----- Calculate PBC Metrics on ${sampleID} -----"
  '''

  { # try

    bedtools bamtobed -bedpe \
    -i !{tmp_bams[0]} \
    | awk 'BEGIN{OFS="\\t"}{print $1,$2,$4,$6,$9,$10}' \
    | grep -v 'MT' | sort | uniq -c \
    | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} \
    ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
    END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > !{sampleID}.pbc.qc
    
  } || { # catch

      echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA" > !{sampleID}.pbc.qc
  
  }


  '''
}
