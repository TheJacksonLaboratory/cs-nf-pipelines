process FRIP_SCORE {
    tag "${ip} vs ${control}"

    cpus 1
    memory 10.GB
    time '10:00:00'


    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? ip+'_vs_'+control+'/macs2' : 'macs2' }", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? ip+'_vs_'+control+'/macs2' : 'macs2' }", pattern: "*.txt", mode: 'copy'

    container 'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0'

    input:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), file(ipbam), val(control), file(controlbam), file(ipflagstat)
    tuple val(ip), file(peak)
    file(peak_count_header)
    file(frip_score_header)

    output:
    tuple val(ip), path("*.tsv"), emit : tsv
    tuple val(ip), path("*.txt"), emit : txt

    script:
    def PEAK_TYPE = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${ip}", \$1 }' | cat $peak_count_header - > ${ip}_peaks.count_mqc.tsv
    READS_IN_PEAKS=\$(intersectBed -a ${ipbam[0]} -b $peak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')i
    grep 'mapped (' $ipflagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${ip}", a/\$1}' | cat $frip_score_header - > ${ip}_peaks.FRiP_mqc.tsv

    find * -type l -name "*.${PEAK_TYPE}" -exec echo -e "macs2/"{}"\\t0,0,178" \\; > ${ip}_peaks.igv.txt

    """
}
