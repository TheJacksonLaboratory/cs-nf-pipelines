process SNPWEIGHTS_INFERANC {
    tag "$sampleID"

    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/snpweights:2.1_bash'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*ancestry.tsv", mode: 'copy'

    input:
    tuple val(sampleID), path(geno), path(ind), path(snp)

    output:
    tuple val(sampleID), path("*ancestry.tsv"), emit: vcf

    script:
    """
    echo "geno: ${geno}" > snpweights.par
    echo "snp: ${snp}" >> snpweights.par
    echo "ind: ${ind}" >> snpweights.par
    echo "snpwt: ${params.snpweights_panel}" >> snpweights.par
    echo "predpcoutput: ${sampleID}.predpc" >> snpweights.par

    inferanc -p snpweights.par 
    printf "nSites\\tAFR\\tEUR\\tEAS\\tAMR\\tSAS\\n" > placeholder.tsv
    cat ${sampleID}.predpc | awk -F' ' '{print \$3,\$10,\$8,\$12,\$9,\$11}' OFS='\\t' >> placeholder.tsv
    cat placeholder.tsv > ${sampleID}.ancestry.tsv
    """
}

// snpweights.par
/*
geno: /path/to/eigenstrat_geno.ext
snp: /path/to/eigenstrat_snp.ext
ind: /path/to/eigenstrat_ind.ext
snpwt: /path/to/snpwt.ext
predpcoutput: maf.predpc
*/
