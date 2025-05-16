// adapted from https://github.com/nf-core/rnaseq/blob/3.14.0/modules/local/rsem_merge_counts/main.nf
process MERGE_RSEM_COUNTS {
    tag "Merge RSEM counts"

    cpus 1  
    memory 24.GB
    time 2.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/py3_perl_pylibs:v2"

    publishDir "${params.pubdir}", pattern: "*rsem.merged.*.tsv", mode:'copy'

    input:
    path("genes/*")
    path("isoforms/*")
    val(prefix)

    output:
    path "*rsem.merged.gene_counts.tsv", emit: gene_counts
    path "*rsem.merged.gene_tpm.tsv", emit: gene_tpm
    path "*rsem.merged.transcript_counts.tsv", emit: transcript_counts
    path "*rsem.merged.transcript_tpm.tsv", emit: transcript_tpm

    script:
    """
    mkdir -p tmp/genes
    cut -f 1,2 `ls ./genes/* | head -n 1` > gene_ids.txt
    for fileid in `ls ./genes/*`; do
        samplename=`basename \$fileid | sed s/\\.genes.results\$//g`
        echo \$samplename > tmp/genes/\${samplename}.counts.txt
        cut -f 5 \${fileid} | tail -n+2 >> tmp/genes/\${samplename}.counts.txt
        echo \$samplename > tmp/genes/\${samplename}.tpm.txt
        cut -f 6 \${fileid} | tail -n+2 >> tmp/genes/\${samplename}.tpm.txt
    done

    mkdir -p tmp/isoforms
    cut -f 1,2 `ls ./isoforms/* | head -n 1` > transcript_ids.txt
    for fileid in `ls ./isoforms/*`; do
        samplename=`basename \$fileid | sed s/\\.isoforms.results\$//g`
        echo \$samplename > tmp/isoforms/\${samplename}.counts.txt
        cut -f 5 \${fileid} | tail -n+2 >> tmp/isoforms/\${samplename}.counts.txt
        echo \$samplename > tmp/isoforms/\${samplename}.tpm.txt
        cut -f 6 \${fileid} | tail -n+2 >> tmp/isoforms/\${samplename}.tpm.txt
    done

    paste gene_ids.txt tmp/genes/*.counts.txt > ${prefix}.rsem.merged.gene_counts.tsv
    paste gene_ids.txt tmp/genes/*.tpm.txt > ${prefix}.rsem.merged.gene_tpm.tsv
    paste transcript_ids.txt tmp/isoforms/*.counts.txt > ${prefix}.rsem.merged.transcript_counts.tsv
    paste transcript_ids.txt tmp/isoforms/*.tpm.txt > ${prefix}.rsem.merged.transcript_tpm.tsv  
    """
}
