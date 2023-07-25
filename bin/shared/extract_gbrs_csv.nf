// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084
def extract_gbrs_csv(csv_file) {

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
            if (!(row.sampleID && row.sex && row.generation && row.fastq_1)){
                log.error "Error in CSV file: Missing field in csv file header. The csv file must have a fields named 'sampleID, sex, generation, fastq_1, {fastq_2}'. Where fastq_2 is optional"
                System.exit(1)
            }
            [row.sampleID.toString(), row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> 

        def meta = [:]

        // Meta data to identify samplesheet
        meta.sampleID = row.sampleID.toString()
        meta.sex = row.sex.toString()
        meta.generation = row.generation.toString()

        // If no lane specified, lane is not considered
        if (row.lane) meta.lane = row.lane.toString()
        else meta.lane = 'NA'

        meta.id = row.sampleID.toString()

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
