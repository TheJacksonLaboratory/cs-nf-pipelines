process ALNTOOLS_BAM2EMASE {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time 10.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
  
    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-547132f'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'alntools' }", pattern: "*.h5", mode: 'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.emase.h5"), emit: emase_h5

    script:
    """
    alntools bam2emase \
                -t ${params.transcripts_info} \
                ${bam} \
                ${bam.baseName}.emase.h5
    """

    stub:
    """
    touch ${bam.baseName}.emase.h5
    """
}

/*
Usage: alntools bam2emase <options> bam_file emase_file

  Convert a BAM file (bam_file) to an EMASE file (emase_file)

Options:
  -c, --chunks INTEGER            number of chunks to process
  -d, --directory DIRECTORY       temp directory
  -m, --mincount INTEGER          minimum count
  --multisample
  -p, --number-processes INTEGER  number of processes
  --rangefile FILE                range file
  -t, --targets FILE              target file
  -v, --verbose                   enables verbose mode
  -h, --help                      Show this message and exit.
*/