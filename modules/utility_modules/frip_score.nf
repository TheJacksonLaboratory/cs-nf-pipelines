process FRIP_SCORE {
    tag "${ip} vs ${control}"

    cpus 1
    memory 10.GB
    time '10:00:00'

    publishDir "${params.pubdir}/${'immuno_precip_samples/' + ip + '_vs_' + control + '/macs2'}", pattern: "*.tsv", mode: 'copy'

    container 'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0'

    input:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), path(ipbam), val(control), path(controlbam), path(ipflagstat), path(peak)
    path(peak_count_header)
    path(frip_score_header)

    output:
    tuple val(ip), path("*.tsv"), emit : tsv

    script:
    def PEAK_TYPE = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${ip}", \$1 }' | cat $peak_count_header - > ${ip}_peaks.count_mqc.tsv
    READS_IN_PEAKS=\$(intersectBed -a ${ipbam[0]} -b $peak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')i
    grep 'mapped (' $ipflagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${ip}", a/\$1}' | cat $frip_score_header - > ${ip}_peaks.FRiP_mqc.tsv
    """
}

/*
IGV steps removed, re-add if IGV is needed: 

    PUBDIR: publishDir "${params.pubdir}/${'comparison/' + ip + '_vs_' + control + '/macs2'}", pattern: "*.txt", mode: 'copy'

    OUTPUT: tuple val(ip), path("*.txt"), emit : txt

    SCRIPT: find * -type l -name "*.${PEAK_TYPE}" -exec echo -e "macs2/"{}"\\t0,0,178" \\; > ${ip}_peaks.igv.txt
*/
