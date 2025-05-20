process CONPAIR_PILEUP {
    tag "$sampleName"

    cpus 1
    memory 20.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/conpair:gatk4.1.5.0_v0.2'

    input:
    tuple val(sampleID), val(sampleName), file(bam), file(bai)
    val(type)

    output:
    tuple val(sampleID), val(sampleName), file("*pileup.txt"), emit: pileup

    script:
    """
    python /Conpair-master/scripts/run_gatk_pileup_for_sample.py -B ${bam} -O ${sampleName}_${type}_pileup.txt --reference ${params.ref_fa} --markers /Conpair-master/conpair/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed
    """
}

// Marker file inputs: `--markers /Conpair-master/data/markers/...` are located in the container. 
