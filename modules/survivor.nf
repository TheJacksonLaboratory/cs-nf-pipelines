process MERGEPACBIO {

    tag "$sampleID"

    cpus 8
    memory { 40.GB * task.attempt }
    time { 4.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/survivor:1.0.7--he513fc3_0'

    input:
        tuple val(sampleID), file(pbsv_vcf)
        tuple val(sampleID), file(sniffles_vcf)
        val(surv_dist)
        val(surv_supp)
        val(surv_type)
        val(surv_strand)
        val(surv_min)
    output:
        tuple val(sampleID), file("${sampleID}_mergedCall.PS.vcf"), emit: survivor_vcf
    script:
        """
        ls ${pbsv_vcf} > vcf_list.txt
        ls ${sniffles_vcf} >> vcf_list.txt
        SURVIVOR merge vcf_list.txt ${surv_dist} ${surv_supp} ${surv_type} ${surv_strand} 0 ${surv_min} ${sampleID}_mergedCall.PS.vcf
        """
}

process MERGEILLUMINA {

    tag "$sampleID"

    cpus 8
    memory { 40.GB * task.attempt }
    time { 4.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/survivor:1.0.7--he513fc3_0'

    input:
        tuple val(sampleID), file(breakdancer_vcf)
        file(delly_vcf)
        file(manta_vcf)
        file(lumpy_vcf)
        val(surv_dist)
        val(surv_supp)
        val(surv_type)
        val(surv_strand)
        val(surv_min)
    output:
        tuple val(sampleID), file("${sampleID}_mergedCall.BDLM.vcf"), emit: survivor_vcf
    script:
        """
        ls ${breakdancer_vcf} > vcf_list.txt
        ls ${delly_vcf} >> vcf_list.txt
        ls ${manta_vcf} >> vcf_list.txt
        ls ${lumpy_vcf} >> vcf_list.txt
        SURVIVOR merge vcf_list.txt ${surv_dist} ${surv_supp} ${surv_type} ${surv_strand} 0 ${surv_min} ${sampleID}_mergedCall.PS.vcf
        """
}

process ANNOTATE{
    tag "$sampleID"

    cpus 1
    memory { 20.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry'    

    publishDir params.pubdir, mode:'copy'

    input:
        tuple val(sampleID), file(survivor_vcf)
        val(seqmode)
    output:
        tuple val(sampleID), file("${sampleID}.merged.overlap.annotated.txt"), emit: annot
    script:
        """
        /usr/bin/env bash ${projectDir}/bin/surv_annot.sh ${sampleID} ${survivor_vcf} ${seqmode}
        """
}

process SUMMARIZE{

    tag "$sampleID"

    cpus 1
    memory { 20.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry' 

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    input:
        tuple val(sampleID), file(survivor_vcf)
    output:
        path("${sampleID}.survivor_summary.csv"), emit: summary
    script:
        """
        /usr/bin/env python ${projectDir}/bin/sv_to_table.py -v ${survivor_vcf} -o ${sampleID}.survivor_summary.csv
        """
}

process PREPBEDS{
    tag "$sampleID"

    cpus 1
    memory { 100.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'rocker/tidyverse:4.2.1'

    input:
        tuple val(sampleID), file(annot)
        file(summary)
    output:
        tuple val(sampleID), file("${sampleID}.ins.bed"), file("${sampleID}.del.bed"), file("${sampleID}.dup.bed"), file("${sampleID}.inv.bed"), file("${sampleID}.tra.bed"), emit: sv_beds
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/surv_annot_process.R ${sampleID}
        """
}

process INTERSECTBEDS{
    tag "$sampleID"

    cpus 1
    memory { 100.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/jaxcompsci/bedtools-sv_refs:2.30.0--hc088bd4_0'

    input:
        tuple val(sampleID), file("${sampleID}.ins.bed"),file("${sampleID}.del.bed"),file("${sampleID}.dup.bed"),file("${sampleID}.inv.bed"),file("${sampleID}.tra.bed")
    output:
        tuple val(sampleID), file("${sampleID}.ins.s.bed"), file("${sampleID}.ins.e.bed"), file("${sampleID}.del.s.bed"), file("${sampleID}.del.e.bed"), file("${sampleID}.inv.e.bed"), file("${sampleID}.tra.e.bed"), file("${sampleID}.dup.e.bed"), file("${sampleID}.ins.genes.bed"), file("${sampleID}.del.genes.bed"), file("${sampleID}.inv.genes.bed"), file("${sampleID}.dup.genes.bed"), file("${sampleID}.tra.genes.bed"), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), emit: intersected_beds
        tuple val(sampleID), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), emit: intersected_exons
    script:
        """
        /usr/bin/env bash ${projectDir}/bin/intersect_beds.sh ${sampleID}
        """
}

process SUMMARIZEINTERSECTIONS{
    
    tag "$sampleID"

    cpus 1
    memory { 100.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'rocker/tidyverse:4.2.1'    

    publishDir params.pubdir, mode:'copy'

    input:
        tuple val(sampleID), file("${sampleID}.ins.bed"),file("${sampleID}.del.bed"),file("${sampleID}.dup.bed"),file("${sampleID}.inv.bed"),file("${sampleID}.tra.bed")
        tuple val(sampleID), file("${sampleID}.ins.s.bed"), file("${sampleID}.ins.e.bed"), file("${sampleID}.del.s.bed"), file("${sampleID}.del.e.bed"), file("${sampleID}.inv.e.bed"), file("${sampleID}.tra.e.bed"), file("${sampleID}.dup.e.bed"), file("${sampleID}.ins.genes.bed"), file("${sampleID}.del.genes.bed"), file("${sampleID}.inv.genes.bed"), file("${sampleID}.dup.genes.bed"), file("${sampleID}.tra.genes.bed"), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed")
        file("${sampleID}.survivor_summary.csv")
        tuple val(sampleID), file("${sampleID}.merged.overlap.annotated.txt")
    output:
        file("${sampleID}.survivor_joined_results.csv")
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/post_filt.R ${sampleID}
        """
}


process ANNOTATEPACBIO{

    tag "$sampleID"

    cpus 1
    memory { 20.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry'  

    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    publishDir params.pubdir, mode:'copy'

    input:
        tuple val(sampleID), file(survivor_vcf)
        tuple val(sampleID), file("ins.exons.bed"), file("del.exons.bed"), file("dup.exons.bed"), file("tra.exons.bed"), file("inv.exons.bed")
    output:
        file("${sampleID}.mergedCalls.InExon.PS.vcf")
    script:
    """
    /usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v ${survivor_vcf} \
        -i ins.exons.bed -d del.exons.bed \
        -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
        -o ${sampleID}.mergedCalls.InExon.PS.vcf
    """
}

process ANNOTATEILLUMINA{

    tag "$sampleID"

    cpus 1
    memory { 20.GB * task.attempt }
    time { 30.m * task.attempt }
    maxRetries 1
    errorStrategy 'retry'  

    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    publishDir params.pubdir, mode:'copy'

    input:
        tuple val(sampleID), file(survivor_vcf)
        tuple val(sampleID), file("ins.exons.bed"), file("del.exons.bed"), file("dup.exons.bed"), file("tra.exons.bed"), file("inv.exons.bed")
    output:
        file("${sampleID}.mergedCalls.InExon.BDLM.vcf")
    script:
    """
    /usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v ${survivor_vcf} \
        -i ins.exons.bed -d del.exons.bed \
        -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
        -o ${sampleID}.mergedCalls.InExon.BDLM.vcf
    """
}