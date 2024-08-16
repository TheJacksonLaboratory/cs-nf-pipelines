#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include {help} from "${projectDir}/bin/help/cnv_array.nf"
include {param_log} from "${projectDir}/bin/log/cnv_array.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_cnv_array_csv.nf"
include {IAAP_CLI} from "${projectDir}/modules/illumina/iaap_cli.nf"
include {BCFTOOLS_GTC2VCF} from "${projectDir}/modules/bcftools/bcftools_gtct2vcf.nf"
include {BCFTOOLS_QUERY_ASCAT} from "${projectDir}/modules/bcftools/bcftools_query_ascat.nf"
include {ASCAT} from "${projectDir}/modules/ascat/ascat_run.nf"
include {ASCAT_ANNOTATION} from "${projectDir}/modules/ascat/ascat_annotation.nf"


// Help if needed
if (params.help) {
    help()
    exit 0
}

// Log parameter info
param_log()
// Parameter validation
if (!params.bpm_file || !params.egt_file) {
    exit 1, "All parameters (bpm_file, egt_file) are required."
}

if (params.csv_input) {
    ch_input = extract_csv(file(params.csv_input, checkIfExists: true))
} else {
    exit 1, "Workflow requires a CSV manifest. See `--help` for information."   
}

GC_file = file(params.gc_file, checkIfExists: true)
RT_file = file(params.rt_file, checkIfExists: true)

// Main workflow
workflow CNV_ARRAY {
    IAAP_CLI(ch_input)
    BCFTOOLS_GTC2VCF(IAAP_CLI.out.gtc)
    BCFTOOLS_QUERY_ASCAT(BCFTOOLS_GTC2VCF.out.gtc2vcf)
    ASCAT(BCFTOOLS_QUERY_ASCAT.out.baf_lrr)
    ASCAT_ANNOTATION(ASCAT.out.seg_ploidy)
}
