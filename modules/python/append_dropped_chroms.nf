process APPEND_DROPPED_CHROMS {
    tag "$strain"

    cpus 2
    memory 5.GB
    time 60.minutes 
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'hdfgroup/h5py:2.7.0'

    publishDir "${params.pubdir}/g2gtools", pattern: "*.gtf", mode:'copy'

    input:
    tuple val(strain), path(vci), path(tbi), path(unmapped), path(gtf)

    output:
    output:
    tuple val(strain), path("*.gtf"), emit: appended_gtf

    when:
    params.append_chromosomes != 'false' || params.append_chromosomes

    script:
    """
    python ${projectDir}/bin/g2gtools/append_dropped_chroms.py -v ${vci} -u ${unmapped} -g ${gtf} -o ${strain}.${params.genome_version}_DroppedChromAppended.gtf
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.gtf
    """
}

/*
    This script appends all full dropped chromosomes back to the GTF file. Entire chroms are dropped due to lack of variants in the SNP or INDEL file. 
    Example: No variants called on chrM, but chrM should be present in the GTF file for downstream EMASE/GBRS use.
             With no variants present in the G2Gtools convert step, the entire chrM would be dropped into the 'unmapped' file.
             If `append_chromosomes` == true, then all fully missing chromosomes will be added back to the GTF in the 
             convert step. The appended annotations in the GTF will be in the source genome coordinates, 
             as no SNPs/InDELs were present. 
*/
