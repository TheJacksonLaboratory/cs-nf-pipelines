import Logos

logo = new Logo()
println '\n'
println logo.show()

// Boolean isNumber(value) {
//     try {
//         new BigInteger(value)
//         return true
//     } catch (Exception exception) {
//         return false
//     }
// }

def isValidInteger(value) {
    value.toString().isInteger()
}

def param_log(){
    if (params.csv_input) {
    log.info """
    GBRS RUN PARAMETER LOG

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
    --csv_input                     ${params.csv_input}
    --download_data                 ${params.download_data}
    --genome_build                  ${params.genome_build}
    --bowtie_index                  ${params.bowtie_index}
    --transcripts_info              ${params.transcripts_info}
    --gbrs_strain_list              ${params.gbrs_strain_list}
    --gene2transcript_csv           ${params.gene2transcript_csv}
    --full_transcript_info          ${params.full_transcript_info}
    --emase_model                   ${params.emase_model}
    --emission_prob_avecs           ${params.emission_prob_avecs}
    --trans_prob_dir                ${params.trans_prob_dir}
    --gene_position_file            ${params.gene_position_file}
    --genotype_grid                 ${params.genotype_grid}
    --founder_hex_colors            ${params.founder_hex_colors}
    --gbrs_expression_threshold     ${params.gbrs_expression_threshold}
    --gbrs_sigma                    ${params.gbrs_sigma}
    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
    } else if (params.concat_lanes) {

    if (!params.sample_generation || !isValidInteger(params.sample_generation) || params.sample_generation > 100) {
        error "ERROR: '--sample_generation' is empty, not an integer, or greater than 100. Input was: \"${params.sample_generation}\"." 
    }

    if (!params.sample_sex || (params.sample_sex != 'M' && params.sample_sex != 'F')) {
        error "ERROR: '--sample_sex' is empty, or not 'M' or 'F'. Input was: \"${params.sample_sex}\"." 
    }

    log.info """
    GBRS RUN PARAMETER LOG

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
    --gene_position_file            ${params.gene_position_file}
    --genotype_grid                 ${params.genotype_grid}
    --founder_hex_colors            ${params.founder_hex_colors}
    --gbrs_expression_threshold     ${params.gbrs_expression_threshold}
    --gbrs_sigma                    ${params.gbrs_sigma}
    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
    } else {

    if (!params.sample_generation || !isValidInteger(params.sample_generation) || params.sample_generation > 100) {
        error "ERROR: '--sample_generation' is empty, not an integer, or greater than 100. Input was: \"${params.sample_generation}\"." 
    }

    if (!params.sample_sex || (params.sample_sex != 'M' && params.sample_sex != 'F')) {
        error "ERROR: '--sample_sex' is empty, or not 'M' or 'F'. Input was: \"${params.sample_sex}\"." 
    }

    log.info """
    GBRS RUN PARAMETER LOG

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
    --gene_position_file            ${params.gene_position_file}
    --genotype_grid                 ${params.genotype_grid}
    --founder_hex_colors            ${params.founder_hex_colors}
    --gbrs_expression_threshold     ${params.gbrs_expression_threshold}
    --gbrs_sigma                    ${params.gbrs_sigma}
    --keep_intermediate             ${params.keep_intermediate}

    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
    }
}

