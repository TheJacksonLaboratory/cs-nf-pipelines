/*
 * Create IGV session file
 */
process IGV {

    container 'quay.io/biocontainers/python:3.8.3'

    publishDir "${params.pubdir}/igv", mode: 'copy'

    input:
    file(fasta)
    file(bigwigs)
    file(peaks)
    file(consensus_peaks)

    output:
    path "*files.txt"  , emit: txt
    path "*.xml"       , emit: xml

    script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
    """
    cat *.txt > igv_files.txt
    ${projectDir}/bin/chipseq/igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'
    """
}
