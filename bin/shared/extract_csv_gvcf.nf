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
            if (!(row.sampleID) | !(row.outputID) | !(row.gvcf) ){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing field in csv file header. The csv file must have fields: 'sampleID', 'outputID', 'gvcf'" + ANSI_RESET)
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

        /* 
            NOTE: Additional metadata parsing could be added here. This function is a minimal implimentation of a csv parser. 
        */

        meta.id = row.sampleID.toString()
        meta.output = row.outputID.toString()
        /* 
            NOTE: Additional ID parsing could be added here. For example a concatenation of patient and sample, if those fields were added to the csv sheet. 
        */
        meta.size = size
        // defines the number of lanes for each sample. 

        // join meta to gvcf
        try {
            file(row.gvcf, checkIfExists: true)
        }
        catch (Exception e) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "The file: " + row.gvcf + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }

        return [meta.id, meta, row.gvcf]

    
    }
}