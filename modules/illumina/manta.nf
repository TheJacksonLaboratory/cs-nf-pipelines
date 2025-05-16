process MANTA {
    tag "$sampleID"

    cpus = 4
    memory 24.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/manta:v1.5.0'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern:"*.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), path("*candidateSmallIndels.vcf.gz"), path("*candidateSmallIndels.vcf.gz.tbi"), emit: manta_smallindel_vcf_tbi
    tuple val(sampleID), path("*diploidSV.vcf.gz"), path("*diploidSV.vcf.gz.tbi"), emit: manta_diploidsv_tbi
    tuple val(sampleID), path("*somaticSV.vcf.gz"), path("*somaticSV.vcf.gz.tbi"), val(meta), val(normal_name), val(tumor_name), val('manta'), emit: manta_somaticsv_tbi
    tuple val(sampleID), path("*candidateSV.vcf.gz"), path("*candidateSV.vcf.gz.tbi"), emit: manta_candidatesv_tbi


    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    # configure manta
    python /usr/local/bin/configManta.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --referenceFasta ${params.ref_fa} \
    --callRegions ${params.callRegions} \
    --runDir ${sampleID}

    # execute manta
    python ${sampleID}/runWorkflow.py -j ${task.cpus} \
    --mode local \
    --memGb ${my_mem}
    
    mv ${sampleID}/results/variants/candidateSmallIndels.vcf.gz ${sampleID}_manta_candidateSmallIndels.vcf.gz
    mv ${sampleID}/results/variants/candidateSmallIndels.vcf.gz.tbi ${sampleID}_manta_candidateSmallIndels.vcf.gz.tbi
    mv ${sampleID}/results/variants/diploidSV.vcf.gz ${sampleID}_manta_diploidSV.vcf.gz
    mv ${sampleID}/results/variants/diploidSV.vcf.gz.tbi ${sampleID}_manta_diploidSV.vcf.gz.tbi
    mv ${sampleID}/results/variants/somaticSV.vcf.gz ${sampleID}_manta_somaticSV.vcf.gz
    mv ${sampleID}/results/variants/somaticSV.vcf.gz.tbi ${sampleID}_manta_somaticSV.vcf.gz.tbi
    mv ${sampleID}/results/variants/candidateSV.vcf.gz ${sampleID}_manta_candidateSV.vcf.gz
    mv ${sampleID}/results/variants/candidateSV.vcf.gz.tbi ${sampleID}_manta_candidateSV.vcf.gz.tbi
    """
}