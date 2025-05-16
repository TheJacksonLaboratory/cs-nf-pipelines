process COSMIC_ANNOTATION {
    tag "$sampleID"

    cpus 1
    memory 1.GB
    time 5.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    input:
    tuple val(sampleID), file(vcf)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
        """
        python \
        ${projectDir}/bin/pta/add_cancer_gene_census.py \
        ${params.cosmic_cgc} \
        ${vcf} \
        ${sampleID}_germline_vep_cosmic_annotated.vcf
        """
}

// cosmic for 'pta' pipeline comes from: 
// curl -H "Authorization: Basic ADD AUTHORIZATION" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/cancer_gene_census.csv
// the above command provides a URL for curl download
// curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v97/cancer_gene_census.csv?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1672931317&Signature=PK8YAGC%2Bh9veZqc7mIZzywkOSf0%3D" --output cancer_gene_census.csv
