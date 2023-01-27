def param_log(){
    if (params.concat_lanes)
        
    log.info """
    ______________________________________________________

        GENERATE MULTIWAY TRANSCRIPTOME PARAMETER LOG

    --comment: ${params.comment}

    Results Published to: ${params.pubdir}
    ______________________________________________________
    --workflow                      ${params.workflow}
    -w                              ${workDir}
    -c                              ${params.config}
    --sample_folder                 ${params.sample_folder}
    --extension                     ${params.extension}
    --pattern                       ${params.pattern}
    --read_type                     ${params.read_type}
    --concat_lanes                  ${params.concat_lanes}
    --concat_sampleID_delim         ${params.concat_sampleID_delim}
    --concat_sampleID_positions     ${params.concat_sampleID_positions}
    --bowtie_index                  ${params.bowtie_index}
    --transcripts_info              ${params.transcripts_info}
    --gbrs_strain_list              ${params.gbrs_strain_list}
    --gene2transcript_list          ${params.gene2transcript_list}
    --full_transcript_list          ${params.full_transcript_list}
    --emase_model                   ${params.emase_model}

    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}
    ______________________________________________________
    """
    else
        log.info """
    ______________________________________________________

        GENERATE MULTIWAY TRANSCRIPTOME PARAMETER LOG

    --comment: ${params.comment}

    Results Published to: ${params.pubdir}
    ______________________________________________________
    --workflow                      ${params.workflow}
    -w                              ${workDir}
    -c                              ${params.config}
    --sample_folder                 ${params.sample_folder}
    --extension                     ${params.extension}
    --pattern                       ${params.pattern}
    --read_type                     ${params.read_type}
    --concat_lanes                  ${params.concat_lanes}
    --concat_sampleID_delim         "N/A"
    --concat_sampleID_positions     "N/A"
    --bowtie_index                  ${params.bowtie_index}
    --transcripts_info              ${params.transcripts_info}
    --gbrs_strain_list              ${params.gbrs_strain_list}
    --gene2transcript_list          ${params.gene2transcript_list}
    --full_transcript_list          ${params.full_transcript_list}
    --emase_model                   ${params.emase_model}

    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}
    ______________________________________________________
    """
}