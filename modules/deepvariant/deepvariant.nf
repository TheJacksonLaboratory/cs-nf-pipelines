process DEEPVARIANT {
    tag "$sampleID"

    cpus 4
    memory 50.GB
    time {10.hour * task.attempt}
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'retry'}
    maxRetries 2

    container 'docker://google/deepvariant:1.9.0'

    input:
    tuple val(sampleID), file(bam), file(bai), val(chrom), val(sex), val(chrom_x), val(chrom_y)

    output:
    tuple val(sampleID), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_channel
    tuple val(sampleID), path("*.gvcf.gz"), path("*.gvcf.gz.tbi"), optional: true, emit: gvcf_channel

    script:
    
    def haploid_contigs = "--haploid_contigs=${chrom_x},${chrom_y}"
    def gvcf_output = params.run_gvcf ? "--output_gvcf ${sampleID}_${chrom}.gvcf.gz" : ""


    """
    if [ ${sex} == "M" ]; then
        run_deepvariant \
        --model_type WGS \
        --ref ${params.ref_fa} \
        ${haploid_contigs} \
        --reads ${bam} \
        --output_vcf ${sampleID}_${chrom}.vcf.gz \
        --sample_name ${sampleID} \
        --num_shards 4 \
        --regions ${chrom} \
        ${gvcf_output} \
        --verbosity 1 
    fi

    if [[ ${sex} == "F" || ${sex} == "U" || ${sex} == "NA" ]]; then 
        run_deepvariant \
        --model_type WGS \
        --ref ${params.ref_fa} \
        --reads ${bam} \
        --output_vcf ${sampleID}_${chrom}.vcf.gz \
        --sample_name ${sampleID} \
        --num_shards 4 \
        --regions ${chrom} \
        ${gvcf_output} \
        --verbosity 1
    fi
    """
}
