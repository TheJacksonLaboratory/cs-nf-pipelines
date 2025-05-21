process MAKE_VCF_LIST {
    tag "$sampleID"
    
    cpus = 1
    time = '00:05:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'ubuntu:20.04'

    input:
    tuple val(sampleID), val(chroms)
    val(order)

    output:
    tuple val(sampleID), file("*.list"), emit: list

    script:
    // Puts Individual Chromosome Files In Order and Then Into List for MergeVCFs
    // convert paths to strings
    string_list = [] 
    for (int i = 0; i < chroms.size(); i++) {
        string_list.add(chroms[i].toString())
        }
    // find matches and put in final list
    sorted=""
    for (int i = 0; i < order.size(); i++) {
        sorted+=(string_list.find{ it.contains('_'+ order[i] + '.vcf')})+"\n"
        }

    """
    echo "$sorted" > ${sampleID}.list
    """
}
