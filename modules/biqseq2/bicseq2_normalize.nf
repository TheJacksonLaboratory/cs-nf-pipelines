process BICSEQ2_NORMALIZE {
    tag "$sampleID"

    cpus = 1
    memory = 8.GB
    time = '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bicseq2:v3'
    // publishDir "${params.pubdir}/${sampleID}", pattern:".txt", mode:'copy'

    input:
    tuple val(sampleID), path(individual_chr_seq_files), val(meta), val(read_ID), val(read_length), val(insert_size)
    val(fasta_file_list)

    output:
    tuple val(sampleID), path("*.norm.bin.txt"), val(meta), val(read_ID), emit: normalized_output

    script:

    // fasta and mappability are set file lists. mappability is set by by read length of sample. 
    // tempSeqs are the .seq files from the prior step. 
    // tempnormpaths are the output bins. 
    // bicseq2config file is a file with just the list of chroms. 
    // sampleID is sampleID
    // out-file is the configuration file used in the next step. 

    // `bicseq2_config_writer` will sort lists by chromosome name, and omit invalid chr names. 
    // Chromosome names in file names must have `chr` in the name. OR the bicseq2config file must be changed to exclude it. 

    fasta_files = fasta_file_list.collect { "$it" }.join(' ')
    
    read_length = read_length.toInteger()

    if( read_length >= 90 && read_length <= 110) {
        mappability_path = params.mappability_directory + '/100'
    } else if( read_length >= 115 && read_length <= 135) {
        mappability_path = params.mappability_directory + '/125'
    } else if( read_length >= 140 && read_length <= 160 ) {
        mappability_path = params.mappability_directory + '/151'
    } else if( read_length >= 240 && read_length <= 260) {
        mappability_path = params.mappability_directory + '/250'
    } else {
        log.info("\nUnsupported read length " + read_length + " in BicSeq2 normalization. This step is about to fail gracefully.\n\n")
        mappability_path = 'error'
    }

    seq_file_list = individual_chr_seq_files.collect { "$it" }.join(' ')


    """
    if [ "${mappability_path}" = "error" ]; then exit 1; fi

    mappability_file_list=`echo ${mappability_path}`

    python3 \
    ${projectDir}/bin/pta/bicseq2_config_writer.py \
    --fa-files ${fasta_files} \
    --mappability-directory ${mappability_path} \
    --temp-seqs ${seq_file_list} \
    --norm-bicseq2-config ${params.bicseq2_chromList} \
    --sample-id ${read_ID} \
    --out-file configuration_file.txt

    rounded_length=`echo ${insert_size} | awk '{print int(\$1+0.5)}'`

    /NBICseq-norm_v0.2.4/NBICseq-norm.pl \
    -l=${read_length} \
    -s=\${rounded_length} \
    -fig=${sampleID}.GCvsRD.pdf \
    -tmp=${sampleID}.tmp \
    configuration_file.txt \
    ${sampleID}.params.out
    """

    stub:
    """
    touch ${read_ID}_chr1.norm.bin.txt
    touch ${read_ID}_chr2.norm.bin.txt
    touch ${read_ID}_chr3.norm.bin.txt
    touch ${read_ID}_chr4.norm.bin.txt
    touch ${read_ID}_chr5.norm.bin.txt
    touch ${read_ID}_chr6.norm.bin.txt
    touch ${read_ID}_chr7.norm.bin.txt
    touch ${read_ID}_chr8.norm.bin.txt
    touch ${read_ID}_chr9.norm.bin.txt
    touch ${read_ID}_chr10.norm.bin.txt
    touch ${read_ID}_chr11.norm.bin.txt
    touch ${read_ID}_chr12.norm.bin.txt
    touch ${read_ID}_chr13.norm.bin.txt
    touch ${read_ID}_chr14.norm.bin.txt
    touch ${read_ID}_chr15.norm.bin.txt
    touch ${read_ID}_chr16.norm.bin.txt
    touch ${read_ID}_chr17.norm.bin.txt
    touch ${read_ID}_chr18.norm.bin.txt
    touch ${read_ID}_chr19.norm.bin.txt
    touch ${read_ID}_chr20.norm.bin.txt
    touch ${read_ID}_chr21.norm.bin.txt
    touch ${read_ID}_chr22.norm.bin.txt
    touch ${read_ID}_chrX.norm.bin.txt
    touch ${read_ID}_chrY.norm.bin.txt
    """
}
// Stub run is to test lower coverage sample data.
