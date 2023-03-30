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
    --sample_generation             ${params.sample_generation}
    --sample_sex                    ${params.sample_sex}
    --bowtie_index                  ${params.bowtie_index}
    --transcripts_info              ${params.transcripts_info}
    --gbrs_strain_list              ${params.gbrs_strain_list}
    --gene2transcript_csv           ${params.gene2transcript_csv}
    --full_transcript_info          ${params.full_transcript_info}
    --emase_model                   ${params.emase_model}
    --emission_prob_avecs           ${params.emission_prob_avecs}
    --trans_prob_dir                ${params.trans_prob_dir}
    --gbrs_expression_threshold     ${params.gbrs_expression_threshold}
    --gbrs_sigma                    ${params.gbrs_sigma}


    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
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
    --sample_generation             ${params.sample_generation}
    --sample_sex                    ${params.sample_sex}
    --bowtie_index                  ${params.bowtie_index}
    --transcripts_info              ${params.transcripts_info}
    --gbrs_strain_list              ${params.gbrs_strain_list}
    --gene2transcript_csv           ${params.gene2transcript_csv}
    --full_transcript_info          ${params.full_transcript_info}
    --emase_model                   ${params.emase_model}
    --emission_prob_avecs           ${params.emission_prob_avecs}
    --trans_prob_dir                ${params.trans_prob_dir}
    --gbrs_expression_threshold     ${params.gbrs_expression_threshold}
    --gbrs_sigma                    ${params.gbrs_sigma}

    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
}