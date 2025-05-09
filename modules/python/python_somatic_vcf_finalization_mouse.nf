process SOMATIC_VCF_FINALIZATION {
    tag "$sampleID"

    cpus 1
    memory 50.GB
    time '3:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*final.*", mode:'copy'
    publishDir "${params.pubdir}/${sampleID}", pattern: "*supplemental.vcf", mode:'copy'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(normal_name), val(tumor_name)
    val(filtered)

    output:
    tuple val(sampleID), file("*final.vcf"), emit: vcf
    tuple val(sampleID), file("*final.txt"), emit: txt
    tuple val(sampleID), file("*supplemental.vcf"), emit: supp_vcf

    script:

    output_suffix = filtered == 'filtered' ?  'filtered' : 'unfiltered'

    """
    python \
    ${projectDir}/bin/pta/annotate_id.py \
    ${vcf} \
    ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_supplemental.vcf

    python \
    ${projectDir}/bin/pta/make_main_vcf.py \
    ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_supplemental.vcf \
    ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.vcf \
    GRCm39

    python \
    ${projectDir}/bin/pta/make_txt.py \
    --vcf ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.vcf \
    --txt ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.txt \
    --vep-version GRCm39 \
    --tumor ${tumor_name} \
    --normal ${normal_name} 
    """
}

    // # NOTE: Differences here relative to the human version of this module: 

    // 1. the `annotate_id.py` script was designed to adjust cosmic IDs and add those to the 'ID' field. 
    //    However, mouse does not have cosmic databases. 
    //    Script is left in place, as GT is also adjusted, 
    //    and the ID field will not be altered as there is no cosmic data.     
    // 2. The `rename_csq_vcf.py` script was designed to rename 1000G CSQ fields in the VCF.
    //    This is not needed for mouse, and the step was removed above. The 'supplemntal' vcf comes from the annotate now.
    // 3. The `make_main_vcf.py` script has provisions for dealing with mouse data. The expected genome logic was updated in the script.
    //    NB: The build ID must match the ID present in the VEP header '##VEP=... assembly="GRCm39"'
    // 4. The `make_txt.py` script was adjusted for GRCm39 and VEP 110 annotations. 
    // 5. The `make_maf.py` script was removed, as HUGO symbol does not apply to mouse, and MAF format is intended for human data
        // python \
        // ${projectDir}/bin/pta/make_maf.py \
        // --vcf ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.vcf \
        // --maf ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.maf \
        // --library WGS \
        // --vep-version GRCm39 \
        // --tumor ${tumor_name} \
        // --normal ${normal_name} \
        // --ensembl-entrez ${params.ensembl_entrez}
    // FROM CONFIG:
    // params.ensembl_entrez='/projects/compsci/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/annotations/GRCm39_ensemblv110_entrez_id_map.csv' // used in somatic vcf finalization.
