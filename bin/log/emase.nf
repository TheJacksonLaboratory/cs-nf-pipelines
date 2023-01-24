def param_log(){
log.info """
______________________________________________________

    GENERATE MULTIWAY TRANSCRIPTOME PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--extension                     ${params.extension}
--pattern                       ${params.pattern}
--read_type                     ${params.read_type}
--concat_lanes                  ${params.concat_lanes}
--bowtie_index                  ${params.bowtie_index}
--transcripts_info              ${params.transcripts_info}
--gbrs_strain_list              ${params.gbrs_strain_list}
--gene2transcript_list          ${params.gene2transcript_list}
--full_transcript_list          ${params.full_transcript_list}
--gbrs_model                    ${params.gbrs_model}

--keep_intermediate               ${params.keep_intermediate}

Project Directory: ${projectDir}
______________________________________________________
"""
}