process GET_READ_LENGTH {
    tag "$sampleID"

    cpus = 1
    time = '00:05:00'

    container 'ubuntu:20.04'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), env(READ_LENGTH), emit: read_length

    script:
    """
    READ_LENGTH=`zcat ${reads[0]} | head -n 400 | awk 'NR%4==2{m=length(\$0)}{print m}' | sort -n | tail -1`
    """
}