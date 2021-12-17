workflow RNASEQ{
  // Step 1: Qual_Stat
  // Step 2: RSEM
  Channel.of( sampleID, fqR1p, fqR2p ).toList().set { sample_fastqs_ch2 }
  // Step 3: Get Read Group Information
  // Step 4a: Picard Alignment Metrics, part 1
  // Step 4b: Picard Alignment Metrics, part 2
  // Step 5: Summary Stats (Human samples only)
  // Step 6a: GATK Coverage Stats (Human samples only)
  // Step 6b: GATK Coverage Stats (Human samples only)
  // Step 7: Move final files to sample directories (human samples)
  // Step 7: Move final files to sample directories (mouse samples)

  workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = ${workflow.complete}
    wfEnd['Duration']     = ${workflow.duration}
    wfEnd['Exit status']  = ${workflow.exitStatus}
    wfEnd['Success']      = ${workflow.success}
      if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = ${workflow.errorMessage}
        wfEnd['.    Report']  = ${workflow.errorReport}
      }
    Summary.show(wfEnd)
  }
}
