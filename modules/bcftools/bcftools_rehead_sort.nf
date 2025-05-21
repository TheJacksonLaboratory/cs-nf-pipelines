process BCFTOOLS_REHEAD_SORT {
    tag "$sampleID"
    
    cpus = 1
    memory = 10.GB
    time '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/unmerged_calls'}", pattern: "${sampleID}_${caller}_sorted.vcf", mode: 'copy'

    container "quay.io/biocontainers/bcftools:1.15--h0ea216a_2"

    input:
        tuple val(sampleID), file(variants)
        val(caller)
        tuple file(fasta), file(fai)

    output:
        tuple val(sampleID), file("${sampleID}_${caller}_sorted.vcf"), emit: vcf_sort

    script:

        if (caller == "lumpy")
            """
            bcftools view -h ${variants} | grep ^## > reheader.txt
            awk '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' ${fai} >> reheader.txt
            bcftools view -h ${variants} | grep -v ^## >> reheader.txt
            printf "${sampleID}_${caller}\n" > sample_head.txt
            bcftools reheader --header reheader.txt \
                --samples sample_head.txt \
                -o ${sampleID}_${caller}Reheader.vcf \
                ${variants}

            bcftools sort ${sampleID}_${caller}Reheader.vcf \
                -O v \
                -o ${sampleID}_${caller}_sorted.vcf
            """

        else if (caller == "delly_sv" | caller == "delly_cnv" | caller == "manta")
            """
            printf "${sampleID}_${caller}\n" > rehead.txt
            bcftools reheader --samples rehead.txt \
                -o ${sampleID}_${caller}Reheader.bcf \
                ${variants}
            bcftools sort ${sampleID}_${caller}Reheader.bcf \
                -O v \
                -o ${sampleID}_${caller}_sorted.vcf
            """
}