def getLibraryId( file ) {
    // file.split("_")[0]
    // The above returned only the first part of the split string. 

    String.join('_', file.split(params.concat_sampleID_delim)[0..params.concat_sampleID_positions-1])
    // NOTE: this will error if concat_sampleID_positions is greater than number of positions. 
    // An Elvis op should be included to catch this case.
    // something like: `xs.length>2 ? xs.tail().join() : xs[1]`
}