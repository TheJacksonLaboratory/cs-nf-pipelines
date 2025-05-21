process CONPAIR {
    tag "$pairName"

    cpus 1
    memory 4.GB
    time '10:00:00'
    container 'quay.io/jaxcompsci/conpair:gatk4.1.5.0_v0.2'
    errorStrategy 'ignore'

    publishDir "${params.pubdir}/${pairName}", pattern:"*.txt", mode:'copy'
    
    input:
    tuple val(sampleID), val(pairName), file(tumor_pileup), file(normal_pileup)

    output:
    tuple val(pairName), file("*_concordance.txt"), emit: concordance
    tuple val(pairName), file("*_contamination.txt"), emit: contamination

    script:
    """
    python /Conpair-master/scripts/verify_concordance.py -T ${tumor_pileup} -N ${normal_pileup} --outfile ${pairName}_concordance.txt -M /Conpair-master/conpair/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt

    python /Conpair-master/scripts/estimate_tumor_normal_contamination.py -T ${tumor_pileup} -N ${normal_pileup} --outfile ${pairName}_contamination.txt -M /Conpair-master/conpair/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt
    """

    stub:
    """
    touch ${pairName}_concordance.txt
    touch ${pairName}_contamination.txt
    """
}
