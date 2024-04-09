process MAKE_GENOME_FILTER {
    tag "$sampleID"
    publishDir "${params.pubdir}/genome_info", mode: 'copy'

    input:
    tuple val(sampleID), path(fai)
    file(blacklist)

    output:
    path('*.bed'), emit: bed
    path('*.sizes'), emit: sizes

    script: 
    fasta="\$(echo ${fai} | sed 's/.fai//g')"
    blacklist_filter = params.blacklist ? "sortBed -i $blacklist -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
    """
    cut -f 1,2 ${fai} > ${fasta}.sizes
    $blacklist_filter > ${fasta}.include_regions.bed
    """
}
