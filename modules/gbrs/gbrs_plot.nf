process GBRS_PLOT  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'

    container 'quay.io/mikewlloyd/gbrs_test:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.pdf", mode: 'copy'

    input:
    tuple val(sampleID), path(interpolated_genoprobs)

    output:
    tuple val(sampleID), file("*gbrs.plotted.genome.pdf"), emit: genotype_pdf


    script:

    """
    cp ${params.base_ref_index_fai} ref.fa.fai
    cp ${params.founder_hex_colors} founder.hexcolor.info

    gbrs plot \
        -i ${interpolated_genoprobs} \
        -o ${sampleID}.gbrs.plotted.genome.pdf \
        -n ${sampleID} 
    """

    stub:
    """
    touch ${sampleID}.gbrs.plotted.genome.pdf
    """
}

/*
usage: gbrs plot [-h] -i GPBFILE [-o OUTFILE] [-n SAMPLE_NAME]
                 [--num-grids GRID_SIZE] [--xt-max XT_MAX] [--xt-size XT_SIZE]
                 [--grid-width WIDTH]

optional arguments:
  -h, --help            show this help message and exit
  -i GPBFILE, --genoprob GPBFILE
  -o OUTFILE, --outfile OUTFILE
  -n SAMPLE_NAME, --sample-name SAMPLE_NAME
  --num-grids GRID_SIZE
  --xt-max XT_MAX
  --xt-size XT_SIZE
  --grid-width WIDTH

  mwl note: num-grids, xt-max, xt-size, and grid-width are used in setting the X tick marks.

    plot arguments: 
    subparser_plot.add_argument('--num-grids', action='store', dest='grid_size', type=int, default=42586)
    subparser_plot.add_argument('--xt-max', action='store', dest='xt_max', type=int, default=4501)
    subparser_plot.add_argument('--xt-size', action='store', dest='xt_size', type=int, default=475)
    subparser_plot.add_argument('--grid-width', action='store', dest='width', type=float, default=0.01)

    plot code:
    ax.set_xticklabels([ '%dM' % xt for xt in np.arange(0, xt_max*grid_size/1000000, xt_size*grid_size/1000000)])
    pyplot.xticks(np.arange(0, xt_max*width, xt_size*width))
*/
