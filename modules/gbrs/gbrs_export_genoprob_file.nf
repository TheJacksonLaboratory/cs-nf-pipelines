process GBRS_EXPORT_GENOPROBS  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'

    container 'quay.io/mikewlloyd/gbrs_test:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.tsv", mode: 'copy'

    input:
    tuple val(sampleID), path(interpolated_genoprobs)

    output:
    tuple val(sampleID), file("*.gbrs.interpolated.genoprobs.tsv"), emit: interpolated_genoprobs_tsv


    script:

    """
    export-genoprob-file \
        -i ${interpolated_genoprobs} \
        -s ${params.gbrs_strain_list} \
        -g ${params.genotype_grid}
    """

    stub:
    """
    touch ${sampleID}.gbrs.interpolated.genoprobs.tsv
    """
}

/*
    usage: export-genoprob-file [-h] -i GPROBFILE -s STRAINS [-g GRIDFILE] [-f]
                                [-c] [-p] [-o OUTFILE] [-v]

    optional arguments:
    -h, --help            show this help message and exit
    -i GPROBFILE, --genoprob-file GPROBFILE
    -s STRAINS, --strains STRAINS
    -g GRIDFILE, --grid-file GRIDFILE
    -f, --in-numpy-format
    -c, --by-chromosomes
    -p, --by-strains
    -o OUTFILE
    -v, --verbose         Toggle DEBUG verbosity


    // From the script, default param setttings: 
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-i",
            "--genoprob-file",
            action="store",
            dest="gprobfile",
            type=str,
            required=True,
        )
        parser.add_argument(
            "-s", "--strains", action="store", dest="strains", type=str, required=True
        )
        parser.add_argument(
            "-g", "--grid-file", action="store", dest="gridfile", type=str, default=None
        )
        parser.add_argument(
            "-f", "--in-numpy-format", action="store_true"  # file format
        )
        parser.add_argument(
            "-c", "--by-chromosomes", action="store_true"  # row dimension
        )
        parser.add_argument(
            "-p", "--by-strains", action="store_false"  # column dimension
        )
        parser.add_argument(
            "-o", action="store", dest="outfile", type=str, default=None
        )
        parser.add_argument(
            "-v", "--verbose", help="Toggle DEBUG verbosity", action="store_true"
        )
        return parser.parse_args()
*/
