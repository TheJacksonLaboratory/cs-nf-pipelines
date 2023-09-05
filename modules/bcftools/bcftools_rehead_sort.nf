process BCFTOOLS_REHEAD_SORT {
    tag "$sampleID"

    cpus = 1
    memory = 10.GB
    time '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/unmerged_calls' : 'unmerged_calls' }", pattern: "${sampleID}_${caller}Sort.vcf", mode: 'copy'

    container "quay.io/biocontainers/bcftools:1.15--h0ea216a_2"

    input:
        tuple val(sampleID), file(variants)
        val(caller)
        tuple file(fasta), file(fai)

    output:
        tuple val(sampleID), file("${sampleID}_${caller}Sort.vcf"), emit: vcf_sort

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
                -o ${sampleID}_${caller}Sort.vcf
            """

        else if (caller == "delly")
            """
            printf "${sampleID}_${caller}\n" > rehead.txt
            bcftools reheader --samples rehead.txt \
                -o ${sampleID}_${caller}Reheader.bcf \
                ${variants}
            bcftools sort ${sampleID}_${caller}Reheader.bcf \
                -O v \
                -o ${sampleID}_${caller}Sort.vcf
            """

        else if (caller == "manta")
            // Note: the Manta variants don't have FORMAT fields so there is no sample name reheader
            """
            bcftools sort ${variants} \
                -O v \
                -o ${sampleID}_${caller}Sort.vcf                
            """	
}