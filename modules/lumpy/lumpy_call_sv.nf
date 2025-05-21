process LUMPY_CALL_SV {
    tag "$sampleID"
    
    cpus = 1
    memory {bam_bwa_lumpy_sort.size() < 40.GB ? 40.GB : 80.GB}
    time {bam_bwa_lumpy_sort.size() < 40.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--refv0.2.0'
    
    input:
        tuple val(sampleID), path(bam_bwa_lumpy_sort), path(bam_bwa_lumpy_sort_bai), path(split_sorted_bam), path(split_sorted_bai), path(dis_sorted_bam), path(dis_sorted_bai)
    
    output:
        tuple val(sampleID), path("${sampleID}_lumpyOut.vcf"), emit: lumpy_vcf

    shell:
        pairend_distro = "pairend_distro.py"
        histo          = sampleID + "_alignBWA_lumpySort.lib1.histo"
        lumpy_vcf      = sampleID + "_lumpyOut.vcf"
        exclude_regions = params.exclude_regions
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
        '''
}
