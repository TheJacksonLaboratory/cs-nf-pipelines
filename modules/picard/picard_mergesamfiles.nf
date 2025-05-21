process PICARD_MERGESAMFILES {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${type + sampleID + '/bam'}"
    }, pattern: "*.bam", mode: 'copy', enabled: params.keep_intermediate


    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("*.bam"), emit: bam

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    prefix = "${sampleID}.mLb.mkD"
    bam_files = bam.findAll { it.toString().endsWith('.bam') }.sort()
    if (bam_files.size() > 1) {
        """
        picard -Xmx${my_mem}G MergeSamFiles \
        ${'INPUT='+bam_files.join(' INPUT=')} \
        OUTPUT=${sampleID}.sorted.bam \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=tmp
        """
    }else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        """
    }
}
