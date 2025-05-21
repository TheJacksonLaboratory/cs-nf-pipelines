process MERGE_SV {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    input:
        tuple val(sampleID), val(normal_name), val(tumor_name), file(manta_vcf), file(manta_vcf_tbi), val(meta_manta), val(manta), file(gripss_vcf), val(gripss_idx), val(meta_gripss), val(gripss)
        val(chrom_list)

    output:
        tuple val(sampleID), file("${sampleID}.manta_gridss_sv.bed"), val(normal_name), val(tumor_name), emit: merged
        tuple val(sampleID), file("${sampleID}.manta_gridss_sv_supplemental.bed"), val(normal_name), val(tumor_name), emit: merged_suppl
        

    script:
        listOfChroms = chrom_list.collect { "$it" }.join(',')

        """
        Rscript ${projectDir}/bin/pta/merge-caller-vcfs.r \
            --vcf=${manta_vcf},${gripss_vcf} \
            --caller=manta,gridss \
            --tumor=${tumor_name} \
            --normal=${normal_name} \
            --build=GRCh38 \
            --slop=300 \
            --allowed_chr=${listOfChroms} \
            --min_sv_length=500 \
            --out_file=${sampleID}.manta_gridss_sv.bed \
            --out_file_supplemental=${sampleID}.manta_gridss_sv_supplemental.bed
        """
}
