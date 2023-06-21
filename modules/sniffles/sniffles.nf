process SNIFFLES {
    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/unmerged_calls' : 'unmerged_calls'}", pattern: "${sampleID}.sniffles_sorted_prefix.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(bam), file(index)
    output:
        tuple val(sampleID), file("${sampleID}.sniffles_sorted_prefix.vcf"), emit: sniffles_vcf
    shell:
        if(params.tandem_repeats)
            """
            sniffles --input ${bam} --vcf ${sampleID}.sniffles_calls.vcf --tandem-repeats ${params.tandem_repeats} --output-rnames -t ${task.cpus}
            
            bash ${projectDir}/bin/clean_sniffles.sh ${sampleID}
            """
        else
            """
            sniffles --input ${bam} --vcf ${sampleID}.sniffles_calls.vcf --output-rnames -t ${task.cpus}

            bash ${projectDir}/bin/clean_sniffles.sh ${sampleID}            
            """
}