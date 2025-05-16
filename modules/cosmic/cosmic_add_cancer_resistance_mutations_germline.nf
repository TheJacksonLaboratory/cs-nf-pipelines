process COSMIC_CANCER_RESISTANCE_MUTATION_GERMLINE {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time 3.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    input:
    tuple val(sampleID), file(vcf)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    """
    python \
    ${projectDir}/bin/pta/add_cancer_resistance_mutations.py \
    ${params.cosmic_cancer_resistance_muts} \
    ${vcf} \
    ${sampleID}_germline_snv_indel_annotated_supplemental.vcf
    """
}

// cosmic for 'pta' pipeline comes from: 
// curl -H "Authorization: Basic ADD AUTHORIZATION" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicResistanceMutations.tsv.gz
// the above command provides a URL for curl download
// curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v97/CosmicResistanceMutations.tsv.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1672933745&Signature=nQ9AFGONT4rDKfM4UZ1cmN4J%2F%2BM%3D" --output CosmicResistanceMutations.tsv.gz
