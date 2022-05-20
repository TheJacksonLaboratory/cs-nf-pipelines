process TRIM_GALORE {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'CONTAINER_TBD'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.fastq", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.fastq.gz_stat"), emit: quality_stats // WHAT STATS DOES IT EMIT? 
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:
  log.info "----- Trim Galore Running on: ${sampleID} -----"

  if ${params.non_directional} {
    directionality = '--non_directional'
  } 

  if (params.read_type == "SE"){
    paired_end = ''
  }

  if (params.read_type == "PE"){
    paired_end = '--paired'
  }


  if (params.workflow == "rrbs"){
    rrbs_flag = '--rrbs'
  }

  """
    trim_galore ${paired_end} ${rrbs_flag} ${directionality} --length ${params.trimLength} -q ${params.qualThreshold}  --stringency ${params.adapOverlap}  -a ${params.adaptorSeq}  --fastqc ${fq_reads}
  """
}





if( params.skip_trimming ){
    ch_trimmed_reads_for_alignment = ch_read_files_trimming
    ch_trim_galore_results_for_multiqc = Channel.from(false)
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( filename.indexOf("_fastqc") > 0 ) "FastQC/$filename"
                else if( filename.indexOf("trimming_report.txt" ) > 0) "logs/$filename"
                else if( !params.save_trimmed && filename == "where_are_my_files.txt" ) filename
                else if( params.save_trimmed && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_read_files_trimming
        file wherearemyfiles from ch_wherearemyfiles_for_trimgalore.collect()

        output:
        set val(name), file('*fq.gz') into ch_trimmed_reads_for_alignment
        file "*trimming_report.txt" into ch_trim_galore_results_for_multiqc
        file "*_fastqc.{zip,html}"
        file "where_are_my_files.txt"

        script:
        def c_r1 = clip_r1 > 0 ? "--clip_r1 $clip_r1" : ''
        def c_r2 = clip_r2 > 0 ? "--clip_r2 $clip_r2" : ''
        def tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 $three_prime_clip_r1" : ''
        def tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 $three_prime_clip_r2" : ''
        def rrbs = params.rrbs ? "--rrbs" : ''
        def cores = 1
        if(task.cpus){
            cores = (task.cpus as int) - 4
            if (params.single_end) cores = (task.cpus as int) - 3
            if (cores < 1) cores = 1
            if (cores > 4) cores = 4
        }
        if( params.single_end ) {
            """
            trim_galore --fastqc --gzip $reads \
              $rrbs $c_r1 $tpc_r1 --cores $cores
            """
        } else {
            """
            trim_galore --fastqc --gzip --paired $reads \
              $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 --cores $cores
            """
        }
    }
}
