process SNIFFLES {
    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"

    container 'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0'

    publishDir "${params.pubdir}/unmerged_calls", pattern: "${sampleID}.sniffles_sorted_prefix.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(bam), file(index)
    output:
        tuple val(sampleID), file("${sampleID}.sniffles_sorted_prefix.vcf"), emit: sniffles_vcf
    shell:
        if(params.tandem_repeats)
            '''
            sniffles --input !{bam} --vcf !{sampleID}.sniffles_calls.vcf --tandem-repeats !{params.tandem_repeats} --output-rnames -t !{task.cpus}
            
            cat !{sampleID}.sniffles_calls.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > !{sampleID}.sniffles_calls_sorted.vcf

            grep "#" !{sampleID}.sniffles_calls_sorted.vcf > !{sampleID}.sniffles_sorted_prefix.vcf
        
            grep -v "#" !{sampleID}.sniffles_calls_sorted.vcf | \
                sed -e 's/None/0/g' >> !{sampleID}.sniffles_sorted_prefix.vcf;  
            '''
        else
            '''
            sniffles --input !{bam} --vcf !{sampleID}.sniffles_calls.vcf --output-rnames -t !{task.cpus}

            cat !{sampleID}.sniffles_calls.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > !{sampleID}.sniffles_calls_sorted.vcf

            grep "#" !{sampleID}.sniffles_calls_sorted.vcf > !{sampleID}.sniffles_sorted_prefix.vcf
        
            grep -v "#" !{sampleID}.sniffles_calls_sorted.vcf | \
                sed -e 's/None/0/g' >> !{sampleID}.sniffles_sorted_prefix.vcf;             
            '''
}