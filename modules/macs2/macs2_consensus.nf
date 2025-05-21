/*
 * Consensus peaks across samples, create boolean filtering file, SAF file for featureCounts
 */
process MACS2_CONSENSUS {
    tag "${antibody}"
    
    cpus 8
    memory 10.GB
    time '10:00:00'

    publishDir "${params.pubdir}/${'consensusCalling_' + antibody + '/macs2'}", pattern: "*_peaks.*", mode: 'copy'

    container 'quay.io/biocontainers/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0'

    input:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), path(peaks) 

    output:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), path('*.bed') , emit: ano
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(''), val(''), path('*.bed') , emit: bed
    tuple val(antibody), path('*.saf') , emit: saf
    tuple val(antibody), path("*.pdf")          , emit: pdf
    tuple val(antibody), path("*.antibody.txt") , emit: txt
    tuple val(antibody), path("*.boolean.txt")  , emit: boolean_txt
    tuple val(antibody), path("*.intersect.txt"), emit: intersect_txt

    when:
    params.macs_gsize && (replicatesExist || multipleGroups) && !params.skip_consensus_peaks

    script: 
    peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    prefix = "${antibody}.consensus_peaks"
    mergecols = params.narrow_peak ? (2..10).join(',') : (2..9).join(',')
    collapsecols = params.narrow_peak ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
    expandparam = params.narrow_peak ? '--is_narrow_peak' : ''
    """
    sort -T '.' -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
        | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

    ${projectDir}/bin/chipseq/macs2_merged_expand.py \\
        ${prefix}.txt \\
        ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${peak_type}","")} \\
        ${prefix}.boolean.txt \\
        --min_replicates $params.min_reps_consensus \\
        $expandparam

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    ${projectDir}/bin/chipseq/plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf

    echo "${prefix}.bed\t${antibody}/${prefix}.bed" > ${prefix}.antibody.txt

    """

}

/*
IGV steps removed, re-add if IGV is needed: 

    OUTPUT: tuple val(antibody), path("*.bed.igv.txt"), emit: igv_txt


    SCRIPT: find * -type f -name "${prefix}.bed" -exec echo -e "macs2/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
*/
