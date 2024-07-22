def help() {
    println '''
Parameter | Default | Description

--idat_folder | /<PATH> | The directory containing IDAT files.
--output_dir | /<PATH> | The directory to store the output files.
--bpm_file | /<PATH> | The path to the BPM file.
--egt_file | /<PATH> | The path to the EGT file.
-w | /<PATH> | The directory for intermediary files and Nextflow processes. This directory can become quite large. Ensure ample storage.
--help | false | Print this help message and exit.
'''
}

