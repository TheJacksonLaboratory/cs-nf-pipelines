    
    
    
    
    
    
    withName: STAR_FOR_SQUID {
        publishDir = [
            path: { "${params.outdir}/star_for_squid" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
        ]
        ext.args = '--twopassMode Basic \
        --chimOutType SeparateSAMold \
        --chimSegmentMin 20 \
        --chimJunctionOverhangMin 12 \
        --alignSJDBoverhangMin 10 \
        --outReadsUnmapped Fastx \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat'
    }