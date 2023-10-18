process RSEM_ALIGNMENT_EXPRESSION {
  tag "$sampleID"

  cpus 12
  memory 90.GB
  time 24.h
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'rsem' }", pattern: "*stats", mode:'copy', enabled: params.rsem_aligner == "bowtie2"
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'rsem' }", pattern: "*results*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'rsem' }", pattern: "*genome.sorted.ba*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'rsem' }", pattern: "*transcript.sorted.ba*", mode:'copy'

  input:
  tuple val(sampleID), path(reads), val(strand_setting), val(read_length)
  val(rsem_ref_path)
  val(rsem_star_prefix)
  val(rsem_ref_prefix)

  output:
  path "*stats"
  path "*results*"
  tuple val(sampleID), path("rsem_aln_*.stats"), emit: rsem_stats
  tuple val(sampleID), path("*.stat/*.cnt"), emit: rsem_cnt
  tuple val(sampleID), path("*genes.results"), emit: rsem_genes
  tuple val(sampleID), path("*isoforms.results"), emit: rsem_isoforms
  tuple val(sampleID), path("*.genome.bam"), emit: bam
  tuple val(sampleID), path("*.transcript.bam"), emit: transcript_bam
  tuple val(sampleID), path("*.genome.sorted.bam"), path("*.genome.sorted.bam.bai"), emit: sorted_genomic_bam
  tuple val(sampleID), path("*.transcript.sorted.bam"), path("*.transcript.sorted.bam.bai"), emit: sorted_transcript_bam
  tuple val(sampleID), path("*final.out"), emit: star_log, optional: true
 
  script:

  if (strand_setting == "reverse_stranded") {
    prob="--forward-prob 0"
  }

  if (strand_setting == "forward_stranded") {
    prob="--forward-prob 1"
  }

  if (strand_setting == "non_stranded") {
    prob="--forward-prob 0.5"
  }

  if (params.read_type == "PE"){
    frag=""
    stype="--paired-end"
    trimmedfq="${reads[0]} ${reads[1]}"
  }
  if (params.read_type == "SE"){
    frag="--fragment-length-mean 280 --fragment-length-sd 50"
    stype=""
    trimmedfq="${reads[0]}"
  }
  if (params.rsem_aligner == "bowtie2"){
    
    rsem_ref_files = file("${rsem_ref_path}/bowtie2/*").collect { "$it" }.join(' ')

    outbam="--output-genome-bam --sort-bam-by-coordinate"
    seed_length="--seed-length ${params.seed_length}"
    sort_command=''
    index_command=''
    intermediate=''
    star_log=''
  }
  if (params.rsem_aligner == "star") {
    outbam="--star-output-genome-bam --sort-bam-by-coordinate"
    seed_length=""
    samtools_mem = (int)(task.memory.giga / task.cpus) 
    // cast to integer rounding down no matter what. If 'round' is used, memory request will exceed limits. 
    sort_command="samtools sort -@ 6 -m 5G -o ${sampleID}.STAR.genome.sorted.bam ${sampleID}.STAR.genome.bam"
    index_command="samtools index ${sampleID}.STAR.genome.sorted.bam"
    intermediate='--keep-intermediate-files'
    star_log="cp ${sampleID}.temp/*.final.out ./${sampleID}.STAR.Log.final.out && rm -r ${sampleID}.temp"

    read_length = read_length.toInteger()

    if( read_length >= 65 && read_length <= 85) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_75/*").collect { "$it" }.join(' ')
    } else if( read_length >= 90 && read_length <= 110 ) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_100/*").collect { "$it" }.join(' ')
    } else if( read_length >= 115 && read_length <= 135 ) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_125/*").collect { "$it" }.join(' ')
    } else if( read_length >= 140 && read_length <= 160 ) {
        rsem_ref_files = file("${rsem_ref_path}/STAR/${rsem_star_prefix}_150/*").collect { "$it" }.join(' ')
    } else {
        log.info("\nUnsupported read length " + read_length + " in RSEM with STAR. RSEM will now fail gracefully.\n\n")
        rsem_ref_files = 'error'
    }

  }

  """
  if [ "${rsem_ref_files}" = "error" ]; then exit 1; fi

  ln -s -f ${rsem_ref_files} . 

  rsem-calculate-expression -p $task.cpus \
  ${prob} \
  ${stype} \
  ${frag} \
  --${params.rsem_aligner} \
  --append-names \
  ${seed_length} \
  ${outbam} \
  ${trimmedfq} \
  ${rsem_ref_prefix} \
  ${sampleID} \
  ${intermediate} \
  --sort-bam-memory-per-thread 5G \
  2> rsem_aln_${sampleID}.stats

  ${star_log}

  ${sort_command}

  ${index_command}
  """
}
