process MAKE_GENOME_FILTER {
    tag "$fai"
    publishDir "${params.pubdir}/genome", mode: 'copy'

    input:
    file(fai)
    file(blacklist)

    output:
    file('*.bed')
    file('*.sizes')

    script: 
    fasta="\$(echo ${fai} | sed 's/.fai//g')"
    blacklist_filter = params.blacklist ? "sortBed -i $blacklist -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
    """
    cut -f 1,2 ${fai} > ${fasta}.sizes
    $blacklist_filter > ${fasta}.include_regions.bed
    """
}

