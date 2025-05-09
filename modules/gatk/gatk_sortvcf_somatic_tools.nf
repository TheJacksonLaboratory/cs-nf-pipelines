process GATK_SORTVCF {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '05:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern:"*_lancet_merged.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), path(list), val(meta), val(normal_name), val(tumor_name), val(tool)

    output:
    tuple val(sampleID), file("*.vcf.gz"), file("*.tbi"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: vcf_tbi

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = list.collect { "-I $it" }.join(' ')

    if (tool == 'lancet_support' && params.gen_org == 'human') {
        chrom_extract = (list =~ /\w+merged_(chr.+)_h.+/)
        tool_name = "lancet_support_"+chrom_extract[0][1]
        tool = chrom_extract[0][1]
        // for final sort merge of lancet confirm, set 'tool_name' to include chrom.
        // set tool to chrom. These steps are required for tuple build as input to final merge. 

    } else if (tool == 'lancet_support' && params.gen_org == 'mouse') {
        chrom_extract = (list =~ /\w+merged_(.+)_h.+/)
        tool_name = "lancet_support_"+chrom_extract[0][1]
        tool = chrom_extract[0][1]
        // for final sort merge of lancet confirm, set 'tool_name' to include chrom.
        // set tool to chrom. These steps are required for tuple build as input to final merge. 

    } else {
        tool_name = tool
    }

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" SortVcf  \
    -SD ${params.ref_fa_dict} \
    ${inputs} \
    -O ${sampleID}_${tool_name}_merged.vcf
        
    bgzip -f -c ${sampleID}_${tool_name}_merged.vcf > ${sampleID}_${tool_name}_merged.vcf.gz
    tabix ${sampleID}_${tool_name}_merged.vcf.gz  
    """
}
