#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/ancestry.nf"
include {param_log} from "${projectDir}/bin/log/ancestry.nf"
include {SAMTOOLS_INDEX} from "${projectDir}/modules/samtools/samtools_index"
include {BCFTOOLS_MPILEUP} from "${projectDir}/modules/bcftools/bcftools_mpileup"
include {BCFTOOLS_CALL} from "${projectDir}/modules/bcftools/bcftools_call"
include {BCFTOOLS_FILTER} from "${projectDir}/modules/bcftools/bcftools_filter"
include {BCFTOOLS_ANNOTATE} from "${projectDir}/modules/bcftools/bcftools_annotate"
include {VCF2EIGENSTRAT} from "${projectDir}/modules/snpweights/snpweights_vcf2eigenstrat"
include {SNPWEIGHTS_INFERANC} from "${projectDir}/modules/snpweights/snpweights_inferanc"

// help if needed
if (params.help){
    help()
    exit 0
}

// main workflow
workflow ANCESTRY_RUN {
    // log params
    param_log()

    // prepare reads channel
    if (params.csv_input) {
        bam_input = extract_csv(file(params.csv_input, checkIfExists: true))
    } else {
        bam_input = Channel.fromFilePairs("${params.sample_folder}/*.bam", checkExists:true, size:1 )
    }
    SAMTOOLS_INDEX(bam_input)
    ANCESTRY(bam_input.join(SAMTOOLS_INDEX.out.bai))
} 

// main workflow
workflow ANCESTRY {
    take:
        bam_bai_input
    main:
        BCFTOOLS_MPILEUP(bam_bai_input)
        BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.vcf)
        BCFTOOLS_FILTER(BCFTOOLS_CALL.out.vcf)
        BCFTOOLS_ANNOTATE(BCFTOOLS_FILTER.out.vcf)
        VCF2EIGENSTRAT(BCFTOOLS_ANNOTATE.out.vcf)
        SNPWEIGHTS_INFERANC(VCF2EIGENSTRAT.out.eigenstrat_files)
} 

def extract_csv(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!(row.sampleID && row.bam)){
                log.error "Error in CSV file: Missing field in csv file header. The csv file must have fields named: 'sampleID, bam"
                System.exit(1)
            }
            [row.sampleID.toString(), row.bam]
        }
}
