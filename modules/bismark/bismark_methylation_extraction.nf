process BISMARK_ALIGNMENT {
  tag "$sampleID"

  cpus 8
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'CONTAINER_TBD'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.sam", mode:'copy' // WHAT IS OUTPUT NEEDS TO BE DETERMINED. 

  input:
  tuple val(sampleID), file(sam)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam // WHAT IS OUTPUT NEEDS TO BE DETERMINED. 

  script:
  log.info "----- Bismark Alignment Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    paired_end = ''
    }
  if (params.read_type == "PE"){
    paired_end = '--paired-end'
    }

  if ${params.non_directional} {
    directionality = '--non_directional'
  } 

  """
  bismark_methylation_extractor ${paired_end}  
  
  // PICK UP HERE TO SET THE PARAMS. MOVE PARAMS TO CONFIG ETC. 
  --output  {in_2}  -multicore ${task.cpus}  {Comprehensive} {overlap} {ignore_forward} {ignore_reverse} {ignore_3prime_forward} {ignore_3prime_reverse}            {BedGraph}    --report  {in_1}
  """ // NOTE: OUTPUT DIR IS CALLED....IS THIS A NEEDED THING? 
}








    <!-- Files:
        Ins:
          1: bismark sam file
          2: outdir
        Outs:
          1: Methylation output
    -->
    
 
    <option name="overlap"                   command_text=""                 value="--no_overlap" />
    <option name="ignore_forward"            command_text="--ignore"         value="0" />
    <option name="ignore_reverse"            command_text="--ignore_r2"      value="0" />
    <option name="ignore_3prime_forward"     command_text="--ignore_3prime"  value="0" />
    <option name="ignore_3prime_reverse"     command_text="--ignore_3prime_r2"   value="0" />
    <option name="BedGraph"                  command_text=""              value="--bedgraph" />
    <option name="Comprehensive"             command_text=""              value="--comprehensive" />


    <command program="force_gd_graph.pl" />

    <command program="">

     

    </command>


</tool>
