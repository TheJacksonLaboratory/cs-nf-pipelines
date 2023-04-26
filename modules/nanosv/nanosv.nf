process NANOSV {
    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"

    container 'quay.io/jaxcompsci/nanosv:1.2.4'

    publishDir "${params.pubdir}/unmerged_calls", pattern: "${sampleID}.nanosv_sorted_prefix.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(bam), file(index)
    output:
        tuple val(sampleID), file("${sampleID}.nanosv_sorted_prefix.vcf"), emit: nanosv_vcf
    shell:
        '''
        NanoSV -o !{sampleID}.nanosv_calls.vcf -t !{task.cpus} -c !{projectDir}/config/nanosv_config.ini !{bam}

        cat !{sampleID}.nanosv_calls.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > !{sampleID}.nanosv_calls_sorted.vcf

        grep "#" !{sampleID}.nanosv_calls_sorted.vcf > !{sampleID}.nanosv_sorted_prefix.vcf
        
        grep -v "#" !{sampleID}.nanosv_calls_sorted.vcf | \
            awk 'BEGIN{FS=OFS="\t"} {sub("^", "NanoSV.", $3)}; 1' | \
            sed -e 's/None/0/g' >> !{sampleID}.nanosv_sorted_prefix.vcf;   
        '''
}