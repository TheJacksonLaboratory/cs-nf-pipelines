import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

if (params.csv_input)
log.info """
EMASE RUN PARAMETER LOG

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
--keep_intermediate             ${params.keep_intermediate}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.concat_lanes)
log.info """
EMASE RUN PARAMETER LOG

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
--genome_build                  ${params.genome_build}
--bowtie_index                  ${params.bowtie_index}
--transcripts_info              ${params.transcripts_info}
--gbrs_strain_list              ${params.gbrs_strain_list}
--gene2transcript_csv           ${params.gene2transcript_csv}
--full_transcript_info          ${params.full_transcript_info}
--emase_model                   ${params.emase_model}
--keep_intermediate             ${params.keep_intermediate}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
log.info """
EMASE RUN PARAMETER LOG

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
--genome_build                  ${params.genome_build}
--bowtie_index                  ${params.bowtie_index}
--transcripts_info              ${params.transcripts_info}
--gbrs_strain_list              ${params.gbrs_strain_list}
--gene2transcript_csv           ${params.gene2transcript_csv}
--full_transcript_info          ${params.full_transcript_info}
--emase_model                   ${params.emase_model}
--keep_intermediate             ${params.keep_intermediate}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}