process GATK_UPDATEVCFSEQUENCEDICTIONARY {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(vcf), path(tbi), val(meta), val(normal_name), val(tumor_name), val(tool)
    val(input_tool)

    output:
    tuple val(sampleID), file("*.vcf.gz"), file("*.tbi"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: vcf_tbi

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    replace_gq_string = input_tool == 'svaba' ? "sed -i 's/GQ,Number=1,Type=Integer/GQ,Number=1,Type=String/g' ${vcf.baseName}.reheaded.vcf && sed -i 's/PL,Number=G,Type=Integer/PL,Number=.,Type=Float/g' ${vcf.baseName}.reheaded.vcf" : ''

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" UpdateVCFSequenceDictionary  \
    --source-dictionary ${params.ref_fa_dict} \
    -V ${vcf} \
    --replace true \
    -O ${vcf.baseName}.reheaded.vcf

    ${replace_gq_string}

    bgzip -f -c ${vcf.baseName}.reheaded.vcf > ${vcf.baseName}.reheaded.vcf.gz
    tabix ${vcf.baseName}.reheaded.vcf.gz
    """
}
