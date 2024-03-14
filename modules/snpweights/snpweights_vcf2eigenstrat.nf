process VCF2EIGENSTRAT {
    tag "$sampleID"

    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v1'

    input:
    tuple val(sampleID), path(vcf)

    output:
    tuple val(sampleID), path("*.geno"), path("*.ind"), path("*.reformatted.snp"), emit: eigenstrat_files

    script:
    """
    python ${projectDir}/bin/ancestry/vcf2eigenstrat.py -o ${sampleID} -v ${vcf} && \
    paste -d ',' <(cat ${sampleID}.snp | sed -e 's/    /,/g' | cut -f 1,2,3,4 -d ',') <(cat ${sampleID}.snp | sed -e 's/    /,/g' | cut -f 1 -d ',' | cut -f 3,4 -d ':' | sed -e 's/:/,/g') | sed -e 's/,/    /g' > ${sampleID}.reformatted.snp
    """
}
// script from: https://github.com/mathii/gdc/blob/master/vcf2eigenstrat.py
// Note: '    ' is used in the SNP file rather than \t. 
