process CONPAIR {
    tag "$pairName"

    cpus 1
    memory 4.GB
    time '10:00:00'
    container 'quay.io/jaxcompsci/conpair:v0.2'
    errorStrategy 'ignore'

    publishDir "${params.pubdir}/${pairName}", pattern:"*.txt", mode:'copy'
    
    input:
    tuple val(sampleID), val(pairName), file(tumor_pileup), file(normal_pileup)

    output:
    tuple val(pairName), file("*_concordance.txt"), emit: concordance
    tuple val(pairName), file("*_contamination.txt"), emit: contamination

    script:
    """
    python2 /Conpair-0.2/scripts/verify_concordance.py -T ${tumor_pileup} -N ${normal_pileup} --outfile ${pairName}_concordance.txt -M /Conpair-0.2/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt

    python2 /Conpair-0.2/scripts/estimate_tumor_normal_contamination.py -T ${tumor_pileup} -N ${normal_pileup} --outfile ${pairName}_contamination.txt -M /Conpair-0.2/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt
    """

    stub:
    """
    touch ${pairName}_concordance.txt
    touch ${pairName}_contamination.txt
    """
}
