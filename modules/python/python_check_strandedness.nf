
process CHECK_STRANDEDNESS {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '1:00:00'
    errorStrategy 'finish'

    container 'quay.io/jaxcompsci/how-are-we-stranded-here:v1.0.1-e6ce74d'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), env(STRAND), emit: strand_setting

    script:
    paired = params.read_type == 'PE' ? "-r2 ${reads[1]}" : ''

    """
    check_strandedness -g ${params.strandedness_gtf} -k ${params.strandedness_ref} -r1 ${reads[0]} ${paired} > ${sampleID}_strandedness.txt 2>&1

    data_type=`grep "Data is likely" ${sampleID}_strandedness.txt`

    if  [[ \$data_type == *RF* ]] ; then
        STRAND='reverse_stranded'
    elif [[ \$data_type == *FR* ]] ; then
        STRAND='forward_stranded'
    elif [[ \$data_type == *unstranded* ]] ; then
        STRAND='non_stranded'
    else
        echo "RNA Seq data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6). Please check your data for low quality and contaminating reads before proceeding."; exit 1;
    fi

    """
}

// Data is likely RF/fr-firststrand
// Data is likely FR/fr-secondstrand
// Data is likely unstranded
