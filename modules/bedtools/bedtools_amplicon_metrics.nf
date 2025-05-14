process TARGET_COVERAGE_METRICS {
    tag "$sampleID"

    cpus 4
    memory 20.GB    
    time = '10:00:00' 
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*coverage_metrics.txt", mode: 'copy'
    
    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0' // note: version difference over other bedtools modules. The 2.23.0 container was failing to parse bed target file. 

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*coverage_metrics.txt"), emit: qc_metrics

    shell:
    '''
    ### total bases (B) that map/align to the on-target (OT) region
    bases_on_target=$(coverageBed -a !{params.target_gatk} -b !{bam} | awk '{if($7>0) total+=$7}END{print total}')

    ### Total length covered by BAM alignment.
    total_bases_covered=$(genomeCoverageBed -ibam !{bam} -bg | awk '{if($4>0) total += ($3-$2)}END{print total}')

    ## Bot / Btot
    awk -v a="$bases_on_target" -v b="$total_bases_covered" 'BEGIN { printf "on_target_percent\\t%s\\n", (a/b)*100 }' </dev/null > !{sampleID}_amplicon_coverage_metrics.txt

    ## Get the depth of 20% of the Average coverage
    perc_mean=$(coverageBed -d -a !{params.target_gatk} -b !{bam} | awk '{if($7>0) total+=1;s+=$7}END{print (s/total)*.2}')

    ### total capture array bases
    total_target_bases=$(awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' !{params.target_gatk})

    ### compute total bases that exceed 20% coverage in capture target region, and calculated coverage uniformity
    coverageBed -d -a !{params.target_gatk} -b !{bam} | awk -v percmean=$perc_mean -v totalbases=$total_target_bases '{if($7>percmean) total+=1;s+=$7}END{print "coverage_uniformity\\t"(total/totalbases)*100}' >> !{sampleID}_amplicon_coverage_metrics.txt

    '''
}


/*
Calculations Per: https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/application-note/primerclip-a-tool-for-trimming-primer-sequences-application-note.pdf?sfvrsn=cf83e107_14
*/