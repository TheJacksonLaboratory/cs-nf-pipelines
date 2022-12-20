process BICSEQ2_NORMALIZE {
    tag "$sampleID"

    cpus = 1
    memory = 8.GB
    time = '05:00:00'
    errorStrategy 'finish'

    container 'quay.io/jaxcompsci/bicseq2:latest'
    // publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'biqseq2' }", pattern:".txt", mode:'copy'

    input:
    tuple val(sampleID), file(individual_chr_seq_files), val(read_ID), val(read_length), val(insert_size)
    val(fasta_file_list)

    output:
    tuple val(sampleID), file("*.norm.bin.txt"), emit: normalized_output

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
    
    if( read_length == '100') {
        mappability_path = params.mappability_directory + '/100'
    } else if( read_length == '125') {
        mappability_path = params.mappability_directory + '/125'
    } else if( read_length == '150' || read_length == '151' ) {
        mappability_path = params.mappability_directory + '/151'
    } else if( read_length == '250') {
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
    ${projectDir}/bin/sv/bicseq2_config_writer.py \
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
    touch ${sampleID}_chr1.norm.bin.txt
    touch ${sampleID}_chr2.norm.bin.txt
    touch ${sampleID}_chr3.norm.bin.txt
    touch ${sampleID}_chr4.norm.bin.txt
    touch ${sampleID}_chr5.norm.bin.txt
    touch ${sampleID}_chr6.norm.bin.txt
    touch ${sampleID}_chr7.norm.bin.txt
    touch ${sampleID}_chr8.norm.bin.txt
    touch ${sampleID}_chr9.norm.bin.txt
    touch ${sampleID}_chr10.norm.bin.txt
    touch ${sampleID}_chr11.norm.bin.txt
    touch ${sampleID}_chr12.norm.bin.txt
    touch ${sampleID}_chr13.norm.bin.txt
    touch ${sampleID}_chr14.norm.bin.txt
    touch ${sampleID}_chr15.norm.bin.txt
    touch ${sampleID}_chr16.norm.bin.txt
    touch ${sampleID}_chr17.norm.bin.txt
    touch ${sampleID}_chr18.norm.bin.txt
    touch ${sampleID}_chr19.norm.bin.txt
    touch ${sampleID}_chr20.norm.bin.txt
    touch ${sampleID}_chr21.norm.bin.txt
    touch ${sampleID}_chr22.norm.bin.txt
    touch ${sampleID}_chrX.norm.bin.txt
    touch ${sampleID}_chrY.norm.bin.txt
    """

}