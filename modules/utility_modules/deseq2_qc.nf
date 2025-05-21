process DESEQ2_QC {
    tag "${antibody}"
    
    cpus 1
    memory 15.GB
    time '10:00:00'

    container 'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0'

    publishDir "${params.pubdir}/${'consensusCalling_' + antibody + '/deseq2'}", mode: 'copy'

    input:
    tuple val(antibody), path(counts)
    file(deseq2_pca_header)
    file(deseq2_clustering_header)

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*.rds"                , optional:true, emit: rds
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors

    script:
    prefix = "${antibody}.consensus_peaks"
    bam_ext = params.read_type == 'SE'  ? '.mLb.clN.sorted.bam' : '.mLb.clN.bam'
    vst = params.deseq2_vst ? '--vst TRUE' : ''
    peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    ${projectDir}/bin/chipseq/deseq2_qc.r \\
        --count_file $counts \\
        --sample_suffix '$bam_ext' \\
        --outdir ./ \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        --id_col 1 --count_col 7 --vst TRUE

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
    sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
    sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
    """
}
