process LUMPY_CALL_SV {
    tag "$sampleID"
    
    cpus = 1
    memory = 40.GB
    time = "10:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'
    
    input:
        tuple val(sampleID), file(bam_bwa_lumpy_sort), file(bam_bwa_lumpy_sort_bai), file(split_sorted_bam), file(split_sorted_bai), file(dis_sorted_bam), file(dis_sorted_bai)
    
    output:
        tuple val(sampleID), file("${sampleID}_lumpySort.vcf"), emit: lumpy_vcf

    shell:
        pairend_distro = "pairend_distro.py"
        histo          = sampleID + "_alignBWA_lumpySort.lib1.histo"
        lumpy_vcf      = sampleID + "_lumpyOut.vcf"
        lumpy_sort_vcf = sampleID + "_lumpySort.vcf"
        exclude_regions = "/ref_data/mm10.gaps.centro_telo.scafold.exclude.bed"
        '''
        RG_ID=$(samtools view -H !{bam_bwa_lumpy_sort} | grep '^@RG' | sed "s/.*ID:\\([^\\t]*\\).*/\\1/g")
        #orig: metrics=$(samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>&1
        samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort} | tail -n+100000 > pre_metrics 2>/dev/null
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
            -pe id:"${RG_ID}",bam_file:!{dis_sorted_bam},histo_file:!{histo},mean:"${mean}",stdev:"${std_dev}",read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
            -sr id:"${RG_ID}",bam_file:!{split_sorted_bam},back_distance:10,weight:1,min_mapping_threshold:20 \
            > !{lumpy_vcf}

        bash !{projectDir}/bin/vcfSort.sh !{lumpy_vcf} !{lumpy_sort_vcf}
        '''
}