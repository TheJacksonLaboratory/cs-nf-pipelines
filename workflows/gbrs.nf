/*
NOTE: THIS PIPELINE WILL CONTAIN ALL STEPS IN EMASE, WITH THE ADDITON OF THE FOLLOWING COMMANDS

    gbrs reconstruct \
        -e ${sampleID}.multiway.genes.tpm \
        -t ${ref_tranprob}/tranprob.DO.${generation}.${sex}.npz \
        -x ${ref_avecs} \
        -g ${ref_GenePosOrdered} \
        -o ${sampleID}

    gbrs quantify \
        -i ${merged} \
        -G ${sampleID}.genotypes.tsv \
        -g ${ref_Gene2Transcripts} \
        -L ${ref_GbrsHybridTargets} \
        -M ${model} \
        --report-alignment-counts \
        -o ${sampleID}

    gbrs interpolate \
        -i ${sampleID}.genoprobs.npz \
        -g ${ref_GenomeGrid} \
        -p ${ref_GenePosOrdered_ver} \
        -o ${sampleID}.gbrs.interpolated.genoprobs.npz

    gbrs plot \
        -i ${sampleID}.gbrs.interpolated.genoprobs.npz \
        -o ${sampleID}.gbrs.plotted.genome.pdf \
        -n ${sampleID}
    """

}

process Export_Genoprob {
    label 'export'
    publishDir path:sample_outdir, mode:'copy', pattern:"*"

    input:
    path(ref_GenomeGrid)
    path(genoprob) from ch_genoprobs

    output:
    path "*" into export_out

    script:
    """
    export-genoprob-file \
        -i ${genoprob} \
        -s A,B,C,D,E,F,G,H \
        -g ${ref_GenomeGrid}
    """
}


*/