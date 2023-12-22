process MANTA_CALL {
    tag "$sampleID"
    
    cpus = 8
    memory {bam.size() < 40.GB ? 40.GB : 80.GB}
    time {bam.size() < 40.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/manta:1.6.0--py27_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/unmerged_calls' : 'unmerged_calls' }", pattern: "${sampleID}_mantaSort.vcf", mode: 'copy'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
        tuple file(fasta), file(fai)
    
    output:
        //tuple val(sampleID), file("${sampleID}_mantaCandidate.vcf"), emit: manta_sv
        tuple val(sampleID), file("${sampleID}_mantaDiploidSV.vcf"), emit: manta_sv

    script:
        """
        rm -rf mantaSVOut
        python /usr/local/bin/configManta.py \
            --runDir mantaSVOut \
            --bam ${bam} \
            --referenceFasta ${fasta}
        python ./mantaSVOut/runWorkflow.py -m local -j ${task.cpus}
        #mv mantaSVOut/results/variants/candidateSV.vcf.gz ./${sampleID}_mantaCandidate.vcf.gz
        mv mantaSVOut/results/variants/diploidSV.vcf.gz ./${sampleID}_mantaDiploidSV.vcf.gz
        gunzip ${sampleID}_mantaDiploidSV.vcf.gz
        """
}
