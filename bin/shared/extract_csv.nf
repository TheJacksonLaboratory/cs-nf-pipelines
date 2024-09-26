// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084

ANSI_RED = "\u001B[31m";
ANSI_RESET = "\u001B[0m";

def extract_csv(csv_file) {
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

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!(row.sampleID) | !(row.fastq_1)){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing field in csv file header. The csv file must have fields: 'sampleID', 'fastq_1', {fastq_2}'. Where fastq_2 is optional." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }
            [row.sampleID.toString(), row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        if (row.sampleID) meta.sampleID = row.sampleID.toString()

        // If no lane specified, lane is not considered
        if (row.lane) meta.lane = row.lane.toString()
        else meta.lane = 'NA'

        // If no replicate specified, replicate is not considered
        if (row.replicate) meta.replicate = row.replicate.toString()
        else meta.replicate = 'NA'
        

        /* 
            NOTE: Additional metadata parsing could be added here. This function is a minimal implimentation of a csv parser. 
        */

        meta.id = row.sampleID.toString()
        /* 
            NOTE: Additional ID parsing could be added here. For example a concatenation of patient and sample, if those fields were added to the csv sheet. 
        */
        meta.size = size
        // defines the number of lanes for each sample. 

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
            
            return [meta.id, meta, row.fastq_1, row.fastq_2]

        } else {
            if (params.read_type == 'PE') {
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "`--read_type` set to 'PE', but only `fastq_1` found in csv manifest. Correct manifest with `fastq_2`, or set `--read_type SE` and restart run." + ANSI_RESET)
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
                System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
                System.exit(1)
            }

            return [meta.id, meta, row.fastq_1]

        }
    }
}

/*
    // Additional check of sample sheet:
    // 1. Each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def sample2patient = [:]

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!sample2patient.containsKey(row.sample.toString())) {
                sample2patient[row.sample.toString()] = row.patient.toString()
            } else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
                log.error('The sample "' + row.sample.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sample.toString()] + '" in the sample sheet.')
                System.exit(1)
            }
        }
*/