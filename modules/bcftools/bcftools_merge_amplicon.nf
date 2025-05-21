process BCFTOOLS_MERGECALLERS {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.vcf", mode:'copy'

    input:
    tuple val(sampleID), path(vcf_haplotypecaller), path(vcf_freebayes)
    val(output_suffix)

    output:
    tuple val(sampleID), path("*.vcf"), emit: vcf

    script:

    """

    awk -F\$'\\t' -v OFS='\\t' '/^[^#]/ { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8";CALLER=GATK_HC",\$9,\$10;next } {print \$0}' ${vcf_haplotypecaller} > hc_temp.txt

    awk -F\$'\\t' -v OFS='\\t' '/^[^#]/ { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8";CALLER=FREEBAYES",\$9,\$10;next } {print \$0}' ${vcf_freebayes} > fb_temp.txt

    sed -i '1 a\\##INFO=<ID=CALLER,Number=1,Type=String,Description="The variant caller used">\\ ' hc_temp.txt

    sed -i '1 a\\##INFO=<ID=CALLER,Number=1,Type=String,Description="The variant caller used">\\ ' fb_temp.txt

    bgzip -c hc_temp.txt > ${vcf_haplotypecaller}.gz
    bgzip -c fb_temp.txt > ${vcf_freebayes}.gz

    echo ${sampleID}_haplotypecaller > hc.txt
    echo ${sampleID}_freebayes > fb.txt

    bcftools reheader ${vcf_haplotypecaller}.gz -s hc.txt -o ${vcf_haplotypecaller.baseName}_rehead.vcf.gz
    bcftools reheader ${vcf_freebayes}.gz -s fb.txt -o ${vcf_freebayes.baseName}_rehead.vcf.gz

    bcftools index ${vcf_haplotypecaller.baseName}_rehead.vcf.gz
    bcftools index ${vcf_freebayes.baseName}_rehead.vcf.gz
    
    bcftools \
    merge \
    --force-samples \
    --no-version \
    --threads ${task.cpus} \
    -o ${sampleID}_${output_suffix}.vcf \
    -i DP:join,AC:join,AF:join,AN:join,ALT_AF:join,DP_HQ:join,CALLER:join \
    ${vcf_haplotypecaller.baseName}_rehead.vcf.gz ${vcf_freebayes.baseName}_rehead.vcf.gz
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