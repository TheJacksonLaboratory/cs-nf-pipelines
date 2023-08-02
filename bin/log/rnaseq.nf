import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

if (params.rsem_aligner != "bowtie2" && params.rsem_aligner != "star") {
  error "'--rsem_aligner': \"${params.rsem_aligner}\" is not valid, supported options are 'bowtie2' or 'star'" 
}

if (params.gen_org != "mouse" && params.gen_org != "human") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse' or 'human'" 
}

if (params.strandedness != null && params.strandedness != "reverse_stranded" && params.strandedness != "forward_stranded" && params.strandedness != "non_stranded") {
  error "'--strandedness': \"${params.strandedness}\" is not valid, supported options are 'reverse_stranded' or 'forward_stranded' or 'non_stranded'" 
}

if (params.pdx && params.rsem_aligner=='bowtie2')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                   ${params.workflow}
--gen_org                    ${params.gen_org}
--genome_build               ${params.genome_build}
--read_type                  ${params.read_type}
--sample_folder              ${params.sample_folder}
--extension                  ${params.extension}
--pattern                    ${params.pattern}
--concat_lanes               ${params.concat_lanes}
--csv_input                  ${params.csv_input}
--download_data              ${params.download_data}
--organize_by                ${params.organize_by}
--pubdir                     ${params.pubdir}
-w                           ${workDir}
--keep_intermediate          ${params.keep_intermediate}
-c                           ${params.config}
--min_pct_hq_reads           ${params.min_pct_hq_reads}
--seed_length                ${params.seed_length}

--pdx                        ${params.pdx}
--xenome_prefix              ${params.xenome_prefix}

--strandedness_ref           ${params.strandedness_ref}
--strandedness_gtf           ${params.strandedness_gtf}
--stradedness                ${params.strandedness}

--rsem_aligner               ${params.rsem_aligner}

Human specific files: 
--rsem_ref_prefix_human      ${params.rsem_ref_prefix_human}
--rsem_ref_files_human       ${params.rsem_ref_files_human}
--picard_dict_human          ${params.picard_dict_human}
--ref_flat_human             ${params.ref_flat_human}
--ribo_intervals_human       ${params.ribo_intervals_human}
--classifier_table           ${params.classifier_table}

Mouse specific files: 
--rsem_ref_prefix_mouse      ${params.rsem_ref_prefix_mouse}
--rsem_ref_files_mouse       ${params.rsem_ref_files_mouse}
--picard_dict_mouse          ${params.picard_dict_mouse}
--ref_flat_mouse             ${params.ref_flat_mouse}
--ribo_intervals_mouse       ${params.ribo_intervals_mouse}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________

"""
else if (params.pdx && params.rsem_aligner=='star')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                   ${params.workflow}
--gen_org                    ${params.gen_org}
--genome_build               ${params.genome_build}
--read_type                  ${params.read_type}
--sample_folder              ${params.sample_folder}
--extension                  ${params.extension}
--pattern                    ${params.pattern}
--concat_lanes               ${params.concat_lanes}
--csv_input                  ${params.csv_input}
--download_data              ${params.download_data}
--organize_by                ${params.organize_by}
--pubdir                     ${params.pubdir}
-w                           ${workDir}
--keep_intermediate          ${params.keep_intermediate}
-c                           ${params.config}
--min_pct_hq_reads           ${params.min_pct_hq_reads}
--seed_length                ${params.seed_length}

--pdx                        ${params.pdx}
--xenome_prefix              ${params.xenome_prefix}

--strandedness_ref           ${params.strandedness_ref}
--strandedness_gtf           ${params.strandedness_gtf}
--stradedness                ${params.strandedness}

--rsem_aligner               ${params.rsem_aligner}

Human specific files: 
--rsem_ref_prefix_human      ${params.rsem_ref_prefix_human}
--rsem_ref_files_human       ${params.rsem_ref_files_human}
--rsem_star_prefix_human     ${params.rsem_star_prefix_human}
--picard_dict_human          ${params.picard_dict_human}
--ref_flat_human             ${params.ref_flat_human}
--ribo_intervals_human       ${params.ribo_intervals_human}
--classifier_table           ${params.classifier_table}

Mouse specific files: 
--rsem_ref_prefix_mouse      ${params.rsem_ref_prefix_mouse}
--rsem_ref_files_mouse       ${params.rsem_ref_files_mouse}
--rsem_star_prefix_mouse     ${params.rsem_star_prefix_mouse}
--picard_dict_mouse          ${params.picard_dict_mouse}
--ref_flat_mouse             ${params.ref_flat_mouse}
--ribo_intervals_mouse       ${params.ribo_intervals_mouse}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________

"""
else if (params.gen_org=='human' && params.rsem_aligner=='bowtie2')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--gen_org              ${params.gen_org}
--genome_build         ${params.genome_build}
--read_type            ${params.read_type}
--sample_folder        ${params.sample_folder}
--extension            ${params.extension}
--pattern              ${params.pattern}
--concat_lanes         ${params.concat_lanes}
--csv_input            ${params.csv_input}
--download_data        ${params.download_data}
--organize_by          ${params.organize_by}
--pubdir               ${params.pubdir}
-w                     ${workDir}
--keep_intermediate    ${params.keep_intermediate}
-c                     ${params.config}
--min_pct_hq_reads     ${params.min_pct_hq_reads}
--hq_pct               ${params.hq_pct}
--strandedness_ref     ${params.strandedness_ref}
--strandedness_gtf     ${params.strandedness_gtf}
--stradedness          ${params.strandedness}
--seed_length          ${params.seed_length}
--rsem_ref_prefix      ${params.rsem_ref_prefix}
--rsem_ref_files       ${params.rsem_ref_files}
--rsem_aligner         ${params.rsem_aligner}
--picard_dict          ${params.picard_dict}
--ref_flat             ${params.ref_flat}
--ribo_intervals       ${params.ribo_intervals}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.gen_org=='human' && params.rsem_aligner=='star')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--gen_org              ${params.gen_org}
--genome_build         ${params.genome_build}
--read_type            ${params.read_type}
--sample_folder        ${params.sample_folder}
--extension            ${params.extension}
--pattern              ${params.pattern}
--concat_lanes         ${params.concat_lanes}
--csv_input            ${params.csv_input}
--download_data        ${params.download_data}
--organize_by          ${params.organize_by}
--pubdir               ${params.pubdir}
-w                     ${workDir}
--keep_intermediate    ${params.keep_intermediate}
-c                     ${params.config}
--min_pct_hq_reads     ${params.min_pct_hq_reads}
--hq_pct               ${params.hq_pct}
--strandedness_ref     ${params.strandedness_ref}
--strandedness_gtf     ${params.strandedness_gtf}
--stradedness          ${params.strandedness}
--seed_length          ${params.seed_length}
--rsem_ref_prefix      ${params.rsem_ref_prefix}
--rsem_ref_files       ${params.rsem_ref_files}
--rsem_aligner         ${params.rsem_aligner}
--rsem_star_prefix     ${params.rsem_star_prefix}
--picard_dict          ${params.picard_dict}
--ref_flat             ${params.ref_flat}
--ribo_intervals       ${params.ribo_intervals}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.gen_org=='mouse' && params.rsem_aligner=='bowtie2')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--extension                     ${params.extension}
--pattern                       ${params.pattern}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
--strandedness_ref              ${params.strandedness_ref}
--strandedness_gtf              ${params.strandedness_gtf}
--stradedness                   ${params.strandedness}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--picard_dict                   ${params.picard_dict}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.gen_org=='mouse' && params.rsem_aligner=='star')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--extension                     ${params.extension}
--pattern                       ${params.pattern}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
--strandedness_ref              ${params.strandedness_ref}
--strandedness_gtf              ${params.strandedness_gtf}
--stradedness                   ${params.strandedness}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--rsem_star_prefix              ${params.rsem_star_prefix}
--picard_dict                   ${params.picard_dict}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else error "Invalid parameters in '--gen_org': ${params.gen_org} and/or in '--rsem_aligner': ${params.rsem_aligner}. Supported options are 'mouse' or 'human' and 'bowtie2' or 'star'."

}
