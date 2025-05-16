process BEDTOOLS_GENOMECOV {
    tag "$sampleID"

    cpus 2
    memory 30.GB 
    time '04:00:00'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples' : 'immuno_precip_samples') : '' 
        "${params.pubdir}/${type + '/' + sampleID + '/bigwig'}"
    }, pattern: "*.txt", mode: 'copy'

    
    container 'quay.io/jaxcompsci/bedtools-sv_refs:2.30.0--hc088bd4_0'
    

    input:
    tuple val(sampleID), path(bam), path(flagstat)

    output:
    tuple val(sampleID), path("*.bedGraph"), emit: bedgraph
    tuple val(sampleID), path("*.txt"), emit: scale_factor


    script:
    pe_fragment = params.read_type == 'SE' ? '' : '-pc'
    extend = (params.read_type == 'SE' && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${sampleID}.scale_factor.txt
        
    bedtools genomecov -ibam ${bam[0]} -bg -scale \$SCALE_FACTOR $pe_fragment $extend | sort -T '.' -k1,1 -k2,2n > ${sampleID}.bedGraph
    """
}
