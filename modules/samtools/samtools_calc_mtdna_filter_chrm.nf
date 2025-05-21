process CALC_MTDNA_FILTER_CHRM {
    tag "$sampleID"

    cpus 4
    memory 20.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*_mtDNA_Content.txt", mode: 'copy'

    input:
    tuple val(sampleID), file(rmdup_bam_file), file(rmdup_bai_file)

    output:
    tuple val(sampleID), file("*.sorted.rmDup.rmChrM.bam"), emit: rmChrM_bam
    tuple val(sampleID), file("*.sorted.rmDup.rmChrM.bam.bai"), emit: rmChrM_bai
    tuple val(sampleID), file("*_mtDNA_Content.txt"), emit: mtdna_log

    shell:
    // Get Mitochondrial and total read counts, calculate %mtDNA and filter Mitochondrial Reads from bam file 

    mt_name = params.gen_org == 'mouse' ?  'MT' : 'chrM'

    '''
    # Get Mitochondrial Read Counts from bam file 
    mtReads=$(samtools idxstats !{rmdup_bam_file} | grep '!{mt_name}' | cut -f 3)
    
    # Get Total Read Counts from bam file
    totalReads=$(samtools idxstats !{rmdup_bam_file} | awk '{SUM += $3} END {print SUM}')

    if [ $mtReads >0 ]
    then
        mtReads=$(echo $mtReads)
    else
        mtReads=$(echo 0)
    fi

    # Calculate %mtDNA
    echo -e 'sampleID\\tPerc mtDNA\\n'!{sampleID}'\\t'$(bc <<< "scale=2;100*$mtReads/$totalReads") >> !{sampleID}_mtDNA_Content.txt

    # Filter Mitochondrial Reads from bam file
    samtools view -@ !{task.cpus} -h !{rmdup_bam_file} \
    | grep -v !{mt_name} \
    | samtools sort -@ !{task.cpus} -O bam \
    -o !{sampleID}.sorted.rmDup.rmChrM.bam \
    && samtools index !{sampleID}.sorted.rmDup.rmChrM.bam
    '''
}
