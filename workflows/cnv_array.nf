#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include {help} from "${projectDir}/bin/help/cnv_array.nf"
include {param_log} from "${projectDir}/bin/log/cnv_array.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_cnv_array_csv.nf"
include {IAAP_CLI} from "${projectDir}/modules/illumina/iaap_cli.nf"
include {BCFTOOLS_GTC2VCF} from "${projectDir}/modules/bcftools/bcftools_gtct2vcf.nf"

// Help if needed
if (params.help) {
    help()
    exit 0
}

// Log parameter info
param_log()
// Parameter validation
if (!params.bpm_file || !params.egt_file) {
    exit 1, "All parameters (idat_folder, bpm_file, egt_file) are required."
}

if (params.csv_input) {
    ch_input = extract_csv(file(params.csv_input, checkIfExists: true))
} else {
    exit 1, "Workflow requires a CSV manifest. See `--help` for information."   
}

// Extract CSV input
ch_input = extract_csv(file(params.csv_input, checkIfExists: true))

// Main workflow
workflow CNV_ARRAY {
    IAAP_CLI(ch_input)
    BCFTOOLS_GTC2VCF(IAAP_CLI.out.gtc)
    BCFTOOLS_GTC2VCF.out.gtc2vcf.view()
}
