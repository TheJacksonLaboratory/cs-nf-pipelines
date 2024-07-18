#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// import modules
// include {help} from "${projectDir}/bin/help/cnv.nf"
// include {param_log} from "${projectDir}/bin/log/cnv.nf"

// Parameter validation
if (!params.idat_folder || !params.output_dir || !params.bpm_file || !params.egt_file) {
    exit 1, "All parameters (idat_folder, output_dir, bpm_file, egt_file) are required."
}

// main workflow
process IAAP_CLI {

    // container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'
    // errorStrategy 'finish'

    input:
    path idat_folder from params.idat_folder
    path output_dir from params.output_dir
    path bpm_file from params.bpm_file
    path egt_file from params.egt_file

    script:
    """
    mkdir -p $output_dir
    chmod a+w $output_dir

    echo "Running IAAP_CLI with BPM file: $bpm_file and EGT file: $egt_file" > $output_dir/iaap_cli.log

    /usr/local/bin/iaap-cli/iaap-cli gencall \
        $bpm_file \
        $egt_file \
        $output_dir \
        --idat-folder $idat_folder \
        --output-gtc >> $output_dir/iaap_cli.log 2>&1
    """
}

workflow {
    IAAP_CLI()
}