process BOWTIE {
    tag "$sampleID"

    cpus 8
    memory 30.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bowtie-samtools:v1'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sampleID), path(fq_read), val(paired_read_num)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam
    tuple val(sampleID), file("*.log"), emit: bowtie_log

    script:
    """
    set -o pipefail
    
    zcat ${fq_read} \\
    | bowtie -p ${task.cpus} -q -a --best --strata --sam -v 3 -x ${params.bowtie_index} - 2> ${sampleID}.bowtie_${paired_read_num}.log \\
    | samtools view -bS - > ${sampleID}_mapped_${paired_read_num}.bam
    """
    // NOTE: This is hard coded to .gz input files. 

    stub:
    """
    touch ${sampleID}.bowtie_${paired_read_num}.log
    touch ${sampleID}_mapped_${paired_read_num}.bam
    """
}

/*

-m 100

OPTIONS USED SUMMARY: 

-a/--all
Report all valid alignments per read or pair (default: off). Validity of alignments is determined by the alignment policy (combined effects of -n, -v, -l, and -e). If more than one valid alignment exists and the --best and --strata options are specified, then only those alignments belonging to the best alignment “stratum” will be reported. Bowtie is designed to be very fast for small -k but bowtie can become significantly slower if -a/--all is specified. If you would like to use Bowtie with -a, consider building an index with a denser suffix-array sample, i.e. specify a smaller -o/--offrate when invoking bowtie-build for the relevant index (see the Performance tuning section for details).

--best
Make Bowtie guarantee that reported singleton alignments are “best” in terms of stratum (i.e. number of mismatches, or mismatches in the seed in the case of -n mode) and in terms of the quality values at the mismatched position(s). Stratum always trumps quality; e.g. a 1-mismatch alignment where the mismatched position has Phred quality 40 is preferred over a 2-mismatch alignment where the mismatched positions both have Phred quality 10. When --best is not specified, Bowtie may report alignments that are sub-optimal in terms of stratum and/or quality (though an effort is made to report the best alignment). --best mode also removes all strand bias. Note that --best does not affect which alignments are considered “valid” by bowtie, only which valid alignments are reported by bowtie. When --best is specified and multiple hits are allowed (via -k or -a), the alignments for a given read are guaranteed to appear in best-to-worst order in bowtie’s output. bowtie is somewhat slower when --best is specified.

--strata
If many valid alignments exist and are reportable (e.g. are not disallowed via the -k option) and they fall into more than one alignment “stratum”, report only those alignments that fall into the best stratum. By default, Bowtie reports all reportable alignments regardless of whether they fall into multiple strata. When --strata is specified, --best must also be specified.

-S/--sam
Print alignments in SAM format. See the SAM output section of the manual for details. To suppress all SAM headers, use --sam-nohead in addition to -S/--sam. To suppress just the @SQ headers (e.g. if the alignment is against a very large number of reference sequences), use --sam-nosq in addition to -S/--sam. bowtie does not write BAM files directly, but SAM output can be converted to BAM on the fly by piping bowtie’s output to samtools view.

-v <int>
Report alignments with at most <int> mismatches. -e and -l options are ignored and quality values have no effect on what alignments are valid. -v is mutually exclusive with -n.

-m <int>
Suppress all alignments for a particular read or pair if more than <int> reportable alignments exist for it. Reportable alignments are those that would be reported given the -n, -v, -l, -e, -k, -a, --best, and --strata options. Default: no limit. Bowtie is designed to be very fast for small -m but bowtie can become significantly slower for larger values of -m. If you would like to use Bowtie for larger values of -k, consider building an index with a denser suffix-array sample, i.e. specify a smaller -o/--offrate when invoking bowtie-build for the relevant index (see the Performance tuning section for details).
*/
