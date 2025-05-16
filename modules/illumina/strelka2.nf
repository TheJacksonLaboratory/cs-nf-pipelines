process STRELKA2 {
    tag "$sampleID"

    cpus = 4
    memory { normal_bam.size() < 60.GB ? 8.GB : 24.GB }
    time { normal_bam.size() < 60.GB ? '08:00:00' : '15:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/strelka2:v2.9.3'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern:"*.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), path(candidateSmallIndels), path(candidateSmallIndels_tbi)

    output:
    tuple val(sampleID), path("*indels.vcf.gz"), path("*indels.vcf.gz.tbi"), val(meta), val(normal_name), val(tumor_name), val('strelka2_indel'), emit: strelka_indel_vcf_tbi
    tuple val(sampleID), path("*snvs.vcf.gz"), path("*snvs.vcf.gz.tbi"), val(meta), val(normal_name), val(tumor_name), val('strelka2_sv'), emit: strelka_snv_vcf_tbi

    script:

    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    # configure strelka
    python /usr/local/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --callRegions ${params.callRegions} \
    --referenceFasta ${params.ref_fa} \
    --indelCandidates ${candidateSmallIndels} \
    --config ${params.strelka_config} \
    --runDir ${sampleID}

    # execute strelka
    python ${sampleID}/runWorkflow.py \
    --mode local \
    --job ${task.cpus} \
    --memGb ${my_mem}

    mv ${sampleID}/results/variants/somatic.snvs.vcf.gz ${sampleID}_strelka_somatic.snvs.vcf.gz
    mv ${sampleID}/results/variants/somatic.snvs.vcf.gz.tbi ${sampleID}_strelka_somatic.snvs.vcf.gz.tbi
    mv ${sampleID}/results/variants/somatic.indels.vcf.gz ${sampleID}_strelka_somatic.indels.vcf.gz
    mv ${sampleID}/results/variants/somatic.indels.vcf.gz.tbi ${sampleID}_strelka_somatic.indels.vcf.gz.tbi
    """
}