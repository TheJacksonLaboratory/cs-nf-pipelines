process GATKv3_5_VARIANTRECALIBRATOR {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk3:3.5-0'
        
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), file(normal_germline_vcf)
    tuple val(sampleID), file(normal_germline_vcf_index)

    output:
    tuple val(sampleID), file("*.*recal.txt"), emit: normal_germline_recal
    tuple val(sampleID), file("*.*tranches.txt"), emit: normal_germline_tranches
    tuple val(sampleID), file("*.*plot.R.txt"), emit: normal_germline_plot_R

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    
    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /usr/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R ${params.ref_fa} \
    -input ${normal_germline_vcf} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_135.b37.vcf \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
    -mode SNP \
    -tranche 99.6 \
    -recalFile ${sampleID}.recal.txt \
    -tranchesFile ${sampleID}.tranches.txt \
    -rscriptFile ${sampleID}.plots.R.txt
    """
}