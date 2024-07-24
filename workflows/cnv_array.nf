#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include {help} from "${projectDir}/bin/help/cnv_array.nf"
include {param_log} from "${projectDir}/bin/log/cnv_array.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_cnv_array_csv.nf"
include {IAAP_CLI} from "${projectDir}/modules/illumina/iaap_cli"
include {BCFTOOLS_GTC2VCF} from "${projectDir}/modules/bcftools_gtc2vcf.nf"

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

// Main workflow
workflow CNV_ARRAY {
    IAAP_CLI(ch_input)
    IAAP_CLI.out.gtc.view()

    // Define paths for BCFTOOLS_GTC2VCF inputs
    bpm_file = file(params.bpm_file)
    csv_file = file(params.csv_file)
    egt_file = file(params.egt_file)
    gtcs_dir = IAAP_CLI.out.gtc
    fasta_file = file(params.fasta_file)
    tsv_file = file(params.tsv_file)

    // Call BCFTOOLS_GTC2VCF process
    BCFTOOLS_GTC2VCF(bpm_file, csv_file, egt_file, gtcs_dir, fasta_file, tsv_file)

    // result.view { files -> 
    //     println "BCF: ${files[0]}"
    //     println "CSI: ${files[1]}"
    //     println "VCF: ${files[2]}"
    //     println "TSV: ${files[3]}"
    // }
    BCFTOOLS_GTC2VCF.out.gtc.view()
}
