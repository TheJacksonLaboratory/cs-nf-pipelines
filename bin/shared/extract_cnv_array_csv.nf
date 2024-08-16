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
            if (!(row.sampleID) || !(row.idat_red) || !(row.idat_green)) {
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "Missing field in csv file header. The csv file must have fields: 'sampleID', 'idat_red', 'idat_green'." + ANSI_RESET)
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

        if (row.idat_red.substring(row.idat_red.lastIndexOf(System.getProperty("file.separator")) + 1).count("_") > 2){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "The file: " + row.idat_red + " containes more than 2 underscores in the name." + ANSI_RESET)
                System.err.println(ANSI_RED + "IDAT files must have only 2 underscores (i.e., xxx_xxx_Grn.idat and xxx_xxx_Red.idat)." + ANSI_RESET)
                System.err.println(ANSI_RED + "GEO (and others) rename files to have more than 2 (i.e., GSMxxx_xxx_xxx_Red.idat). File names must be adjusted prior to running." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
        }
        if (row.idat_green.substring(row.idat_green.lastIndexOf(System.getProperty("file.separator")) + 1).count("_") > 2){
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "The file: " + row.idat_green + " containes more than 2 underscores in the name." + ANSI_RESET)
                System.err.println(ANSI_RED + "IDAT files must have only 2 underscores (i.e., xxx_xxx_Grn.idat and xxx_xxx_Red.idat)." + ANSI_RESET)
                System.err.println(ANSI_RED + "GEO (and others) rename files to have more than 2 (i.e., GSMxxx_xxx_xxx_Red.idat). File names must be adjusted prior to running." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
        }

        // Metadata to identify samplesheet
        def meta = [:]

        if (row.sampleID) meta.sampleID = row.sampleID.toString()

        if (row.gender != "XY" && row.gender != "XX" && row.gender != ""){
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "Gender must be 'XX', 'XY' or empty. " + row.gender + " was provided, and isn't valid." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "-----------------------------------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }

        if (row.gender == "") {
            meta.gender = 'NA'
        } else {
            meta.gender = row.gender.toString()
        }
    
        // join meta to idat, and check file existence. 
        try {
            file(row.idat_red, checkIfExists: true)
        }
        catch (Exception e) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "The file: " + row.idat_red + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }
        try {
            file(row.idat_green, checkIfExists: true)
        }
        catch (Exception e) {
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.err.println(ANSI_RED + "The file: " + row.idat_green + " does not exist. Use absolute paths, and check for correctness." + ANSI_RESET)
            System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
            System.err.println(ANSI_RED + "---------------------------------------------" + ANSI_RESET)
            System.exit(1)
        }

        return [meta.sampleID, meta, row.idat_red, row.idat_green]

    }
}
