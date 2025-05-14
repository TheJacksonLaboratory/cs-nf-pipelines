process BEDTOOLS_SUBTRACT {
    tag "$sampleID"

    cpus = 1
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    publishDir "${params.pubdir}/${sampleID + '/sequenza_cnv'}", pattern:"*segments_naWindowFiltered.txt", mode:'copy'

    input:
    tuple val(sampleID), path(segments), path(na_win)

    output:
    tuple val(sampleID), val('to_maintain_cardinality_for_next_module'), file("*segments_naWindowFiltered.txt"), emit: segments_filtered

    script:

    """
    bedtools subtract -N -f 0.5 -a ${segments} -b ${na_win} > ${sampleID}_tmp_sub.txt

    cat ${projectDir}/bin/wes/sequenza_header.txt ${sampleID}_tmp_sub.txt > ${sampleID}_segments_naWindowFiltered.txt
    """
}

/*
	-A	Remove entire feature if any overlap.  That is, by default,
		only subtract the portion of A that overlaps B. Here, if
		any overlap is found (or -f amount), the entire feature is removed.

	-N	Same as -A except when used with -f, the amount is the sum
		of all features (not any single feature).

    -f	Minimum overlap required as a fraction of A.
        - Default is 1E-9 (i.e., 1bp).
        - FLOAT (e.g. 0.50)

    Note: header file is constant of Sequenza. If header is ever altered by sequenza devs, the header file will require an update.
*/
