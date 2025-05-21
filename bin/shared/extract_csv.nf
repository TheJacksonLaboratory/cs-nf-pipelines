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
            System.err.println(ANSI_RED + "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, and at least one sample." + ANSI_RESET)
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }

        // Reopen the file to read and check the header
        def headerLine
        file(csv_file).withReader('UTF-8') { headerReader ->
            headerLine = headerReader.readLine()
        }
        def headers = headerLine.split(',').collect { it.trim() }
        def requiredHeaders = ['sampleID', 'fastq_1']

        if (params.read_type == 'PE') {
            requiredHeaders << 'fastq_2'
        }
        if (params.merge_inds) {
            requiredHeaders << 'ind'
        }
        if (params.deepvariant) {
            requiredHeaders << 'sex'
        }
        if (params.merge_replicates) {
            requiredHeaders << 'replicate'
        }

        def requiredHeadersStr = requiredHeaders.collect { "'${it}'" }.join(', ')
  
        def missingHeaders = requiredHeaders.findAll { !headers.contains(it) }
        if (missingHeaders) {
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "Missing required header(s) in CSV file: ${missingHeaders.join(', ')}" + ANSI_RESET)
            System.err.println(ANSI_RED + "The csv file must have fields: ${requiredHeadersStr}" + ANSI_RESET)
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }
    }

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!(row.sampleID) | !(row.fastq_1)){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing row data in field: 'sampleID' or 'fastq_1'. These fields can not be empty." + ANSI_RESET)
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
        .map{ row, numLanes ->

        def meta = [:]

        // Meta data to identify samplesheet
        if (row.sampleID) meta.sampleID = row.sampleID.toString()

        // If no lane specified, lane is not considered
        if (row.lane) meta.lane = row.lane.toString()
        else meta.lane = 'NA'

        // If no replicate specified, replicate is not considered
        if (row.replicate) meta.replicate = row.replicate.toString()
        else meta.replicate = 'NA'

        // Parse optional "ind" field
        if (row.ind) meta.ind = row.ind.toString()
        else meta.ind = 'NA'

        // Parse optional "sex" field
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'
        
        // Define the ID field for the sample.
        meta.id = row.sampleID.toString()
        
        // defines the number of lanes for each sample. 
        meta.size = size

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
