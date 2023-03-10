process LUMPY_CALL_SV {
    tag "$sample_name"
    
    cpus = 1
    memory = 40.GB
    time = "10:00:00"

    container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'
    
    input:
        tuple val(sampleID), file(bam_bwa_lumpy_sort), file(split_sorted_bam), file(dis_sorted_bam)
    
    output:
        tuple val(sampleID), file("${sampleID}_lumpySort.vcf"), emit: lumpy_vcf

    shell:
        pairend_distro = "pairend_distro.py"
        histo          = sample_name + "_alignBWA_lumpySort.lib1.histo"
        lumpy_vcf      = sample_name + "_lumpyOut.vcf"
        lumpy_sort_vcf = sample_name + "_lumpySort.vcf"
        exclude_regions = "/ref_data/mm10.gaps.centro_telo.scafold.exclude.bed"
        '''
        RG_ID=$(samtools view -H !{bam_bwa_lumpy_sort[1]} | grep '^@RG' | sed "s/.*ID:\\([^\\t]*\\).*/\\1/g")
        #orig: metrics=$(samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>&1
        samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 > pre_metrics 2>/dev/null
        metrics=$(cat pre_metrics | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>/dev/null \
            && [ $? = 141 ] && echo 'metrics to pairend_distro had exitcode: '$?;
        mean=$(echo "${metrics}" | cut -d " " -f 1)
        mean=$(echo "${mean}"    | cut -d ":" -f 2)
        std_dev=$(echo "${metrics}" | cut -d " " -f 2)
        std_dev=$(echo "${std_dev}" | cut -d ":" -f 2)
        rm pre_metrics;

        lumpy \
            -mw 4 \
            -x !{exclude_regions} \
            -pe id:"${RG_ID}",bam_file:!{dis_sorted_bam[1]},histo_file:!{histo},mean:"${mean}",stdev:"${std_dev}",read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
            -sr id:"${RG_ID}",bam_file:!{split_sorted_bam[1]},back_distance:10,weight:1,min_mapping_threshold:20 \
            > !{lumpy_vcf}

        vcfSort.sh !{lumpy_vcf} !{lumpy_sort_vcf}
        '''
}