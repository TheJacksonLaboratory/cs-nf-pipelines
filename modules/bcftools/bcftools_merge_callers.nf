process BCFTOOLS_MERGECALLERS {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), file(vcf), file(idx), val(meta), val(chrom)

    output:
    tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

    script:

    """
    bcftools \
    merge \
    -r ${chrom} \
    --force-samples \
    --no-version \
    --threads ${task.cpus} \
    -f PASS,SUPPORT \
    -F x \
    -m none \
    -o ${sampleID}_mergedCallers_${chrom}.vcf \
    -i called_by:join,num_callers:sum,MNV_ID:join,supported_by:join \
    ${vcf}
    """
}


/*

About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.
         Note that only records from different files can be merged, never from the same file. For
         "vertical" merge take a look at "bcftools norm" instead.
Usage:   bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]

Options:
        --force-samples                resolve duplicate sample names
        --print-header                 print only the merged header and exit
        --use-header <file>            use the provided header
    -0  --missing-to-ref               assume genotypes at missing sites are 0/0
    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. "PASS,.")
    -F, --filter-logic <x|+>           remove filters if some input is PASS ("x"), or apply all filters ("+") [+]
    -g, --gvcf <-|ref.fa>              merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max
    -i, --info-rules <tag:method,..>   rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [DP:sum,DP4:sum]
    -l, --file-list <file>             read file names from the file
    -m, --merge <string>               allow multiallelic records for <snps|indels|both|all|none|id>, see man page for details [both]
        --no-version                   do not append version and command line to the header
    -o, --output <file>                write output to a file [standard output]
    -O, --output-type <b|u|z|v>        'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
    -r, --regions <region>             restrict to comma-separated list of regions
    -R, --regions-file <file>          restrict to regions listed in a file
        --threads <int>                number of extra output compression threads [0]

*/