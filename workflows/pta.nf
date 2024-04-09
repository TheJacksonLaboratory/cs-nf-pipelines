#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/pta.nf"
include {param_log} from "${projectDir}/bin/log/pta.nf"

include {HS_PTA} from "${projectDir}/subworkflows/hs_pta"
include {MM_PTA} from "${projectDir}/subworkflows/mm_pta"

include {CONCATENATE_PTA_FASTQ} from "${projectDir}/subworkflows/concatenate_pta_fastq"


// help if needed
if (params.help){
    help()
    exit 0
}

// log paramiter info
param_log()

// main workflow
workflow PTA {

    if (params.csv_input) {
        ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
        // Concat local Fastq files from CSV input if required.
            CONCATENATE_PTA_FASTQ(ch_input_sample)
            
    }

    if (params.gen_org == "mouse") {
        
        MM_PTA(CONCATENATE_PTA_FASTQ.out.read_meta_ch)

    } else if (params.gen_org == "human") {

        HS_PTA(CONCATENATE_PTA_FASTQ.out.read_meta_ch)

    } else {
        
    }

}

// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084
def extract_csv(csv_file) {
    ANSI_RED = "\u001B[31m";
    ANSI_RESET = "\u001B[0m";

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines." + ANSI_RESET)
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }
    }

    // Additional check of sample sheet:
    // 1. Each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def patient_sample_lane_combinations_in_samplesheet = []
    def sample2patient = [:]

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!sample2patient.containsKey(row.sampleID.toString())) {
                sample2patient[row.sampleID.toString()] = row.patient.toString()
            } else if (sample2patient[row.sampleID.toString()] != row.patient.toString()) {
                System.err.println(ANSI_RED + 'The sample "' + row.sampleID.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sampleID.toString()] + '" in the sample sheet.' + ANSI_RESET)
                System.exit(1)
            }
        }

    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0

    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all++
            if (!(row.patient && row.sampleID)){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing field in csv file header. The csv file must have a fields: 'patient', 'sampleID', 'lane', 'fastq_1', 'fastq_2'." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            [[row.patient.toString(), row.sampleID.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sampleID)  meta.sampleID  = row.sampleID.toString()

        // If no sex specified, sex is not considered
        // sex is only mandatory for somatic CNV
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 0) sample_count_normal++
        else sample_count_tumor++

        // join meta to fastq
        if (row.fastq_2) {
            if (params.read_type == 'SE') {
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "fastq_2 found in CSV manifest, but `--read_type` set to 'SE'. Set `--read_type PE` and restart run." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            try {
                file(row.fastq_1, checkIfExists: true)
            }
            catch (Exception e) {
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "The file: " + row.fastq_1 + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            try {
                file(row.fastq_2, checkIfExists: true)
            }
            catch (Exception e) {
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "The file: " + row.fastq_2 + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }

            meta.id         = "${row.patient}--${row.sampleID}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)

            meta.size       = 1 // default number of splitted fastq

            return [meta.id, meta, [fastq_1, fastq_2]]

        } else {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "Missing or unknown field in csv file header. Please check your samplesheet" + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }
    }
}