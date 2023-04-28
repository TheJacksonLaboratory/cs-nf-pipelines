// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084
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
            if (!(row.sampleID)){
                log.error "Missing field in csv file header. The csv file must have a field named 'sampleID'."
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

            return [meta.id, meta, row.fastq_1, row.fastq_2]

        } else {
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