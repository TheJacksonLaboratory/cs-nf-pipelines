process JVARKIT_COVERAGE_CAP {
    tag "$sampleID"

    cpus = 1
    memory = 90.GB
    time = '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/jvarkit_samtools:v1'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    jvarkit -Xmx${my_mem}G sortsamrefname --regions ${params.primary_chrom_bed} --samoutputformat BAM ${bam}  |\
    jvarkit -Xmx${my_mem}G biostar154220 --regions ${params.primary_chrom_bed} -d ${params.coverage_cap} --samoutputformat BAM |\
    samtools sort -T tmp -o ${sampleID}_coverageCap.bam -
    """
}

/*
Usage: java -jar dist/jvarkit.jar biostar154220  [options] Files

Usage: biostar154220 [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -d, -n, --depth
      expected coverage.
      Default: 20
    -filter, --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --keep-unmapped
      write unmapped reads
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --query-sorted
      Input was sorted on query name but I promess there is one and only one 
      chromosome: e.g: samtools view -h in.bam 'chr1:234-567' | samtools sort 
      -n -) .
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit
*/
