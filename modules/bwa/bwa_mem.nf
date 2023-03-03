process BWA_MEM {
    tag "${sampleID}"

    cpus 8
    memory 250.GB
    time '72:00:00'
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_6'
    stageInMode 'copy'

    input:
        tuple val(sampleID), file(fq_reads), file(read_group)
        file(bwa_reference)

    output:
        tuple val(sampleID), file(${sampleID}.sam), emit: sam

    script:
    
        if (!params.fastq2)  {
            inputfq="${fq_reads[0]}"
        }
        else {
            inputfq="${fq_reads[0]} ${fq_reads[1]}"
        }
        """
        bwa mem -K 100000000 -R \$(cat $read_group) -t ${task.cpus} -M ${bwa_referece} ${inputfq} > ${sample_name}.sam

        """