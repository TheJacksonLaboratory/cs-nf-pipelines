 */
if( params.aligner =~ /bismark/ ){
    process bismark_align {
        tag "$name"
        publishDir "${params.outdir}/bismark_alignments", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( filename.indexOf(".fq.gz") > 0 ) "unmapped/$filename"
                else if( filename.indexOf("report.txt") > 0 ) "logs/$filename"
                else if( (!params.save_align_intermeds && !params.skip_deduplication && !params.rrbs).every() && filename == "where_are_my_files.txt" ) filename
                else if( (params.save_align_intermeds || params.skip_deduplication || params.rrbs).any() && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_trimmed_reads_for_alignment
        file index from ch_bismark_index_for_bismark_align.collect()
        file wherearemyfiles from ch_wherearemyfiles_for_bismark_align.collect()
        file knownsplices from ch_splicesites_for_bismark_hisat_align.collect().ifEmpty([])

        output:
        set val(name), file("*.bam") into ch_bam_for_bismark_deduplicate, ch_bam_for_bismark_summary, ch_bam_for_preseq
        set val(name), file("*report.txt") into ch_bismark_align_log_for_bismark_report, ch_bismark_align_log_for_bismark_summary, ch_bismark_align_log_for_multiqc
        file "*.fq.gz" optional true
        file "where_are_my_files.txt"

        script:
        // Paired-end or single end input files
        input = params.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

        // Choice of read aligner
        aligner = params.aligner == "bismark_hisat" ? "--hisat2" : "--bowtie2"

        // Optional extra bismark parameters
        splicesites = params.aligner == "bismark_hisat" && params.known_splices ? "--known-splicesite-infile <(hisat2_extract_splice_sites.py ${knownsplices})" : ''
        pbat = params.pbat ? "--pbat" : ''
        non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
        unmapped = params.unmapped ? "--unmapped" : ''
        mismatches = params.relax_mismatches ? "--score_min L,0,-${params.num_mismatches}" : ''
        soft_clipping = params.local_alignment ? "--local" : ''
        minins = bismark_minins ? "--minins $bismark_minins" : ''
        maxins = bismark_maxins ? "--maxins $bismark_maxins" : ''

        // Try to assign sensible bismark memory units according to what the task was given
        multicore = ''
        if( task.cpus ){
            // Numbers based on recommendation by Felix for a typical mouse genome
            if( params.single_cell || params.zymo || params.non_directional ){
                cpu_per_multicore = 5
                mem_per_multicore = (18.GB).toBytes()
            } else {
                cpu_per_multicore = 3
                mem_per_multicore = (13.GB).toBytes()
            }
            // Check if the user has specified this and overwrite if so
            if(params.bismark_align_cpu_per_multicore) {
                cpu_per_multicore = (params.bismark_align_cpu_per_multicore as int)
            }
            if(params.bismark_align_mem_per_multicore) {
                mem_per_multicore = (params.bismark_align_mem_per_multicore as nextflow.util.MemoryUnit).toBytes()
            }
            // How many multicore splits can we afford with the cpus we have?
            ccore = ((task.cpus as int) / cpu_per_multicore) as int
            // Check that we have enough memory, assuming 13GB memory per instance (typical for mouse alignment)
            try {
                tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
                mcore = (tmem / mem_per_multicore) as int
                ccore = Math.min(ccore, mcore)
            } catch (all) {
                log.debug "Warning: Not able to define bismark align multicore based on available memory"
            }
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }

        // Main command
        """
        bismark $input \\
            $aligner \\
            --bam $pbat $non_directional $unmapped $mismatches $multicore $minins $maxins \\
            --genome $index \\
            $soft_clipping \\
            $splicesites
        """
    }

    /*
     * STEP 4 - Bismark deduplicate
     */
    if( params.skip_deduplication || params.rrbs ) {
        ch_bam_for_bismark_deduplicate.into { ch_bam_dedup_for_bismark_methXtract; ch_bam_dedup_for_qualimap }
        ch_bismark_dedup_log_for_bismark_report = Channel.from(false)
        ch_bismark_dedup_log_for_bismark_summary = Channel.from(false)
        ch_bismark_dedup_log_for_multiqc  = Channel.from(false)
    } else {
        process bismark_deduplicate {
            tag "$name"
            publishDir "${params.outdir}/bismark_deduplicated", mode: params.publish_dir_mode,
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            set val(name), file(bam) from ch_bam_for_bismark_deduplicate

            output:
            set val(name), file("*.deduplicated.bam") into ch_bam_dedup_for_bismark_methXtract, ch_bam_dedup_for_qualimap
            set val(name), file("*.deduplication_report.txt") into ch_bismark_dedup_log_for_bismark_report, ch_bismark_dedup_log_for_bismark_summary, ch_bismark_dedup_log_for_multiqc

            script:
            fq_type = params.single_end ? '-s' : '-p'
            """
            deduplicate_bismark $fq_type --bam $bam
            """
        }
    }

    /*
     * STEP 5 - Bismark methylation extraction
     */
    process bismark_methXtract {
        tag "$name"
        publishDir "${params.outdir}/bismark_methylation_calls", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( filename.indexOf("splitting_report.txt" ) > 0 ) "logs/$filename"
                else if( filename.indexOf("M-bias" ) > 0) "m-bias/$filename"
                else if( filename.indexOf(".cov" ) > 0 ) "methylation_coverage/$filename"
                else if( filename.indexOf("bedGraph" ) > 0 ) "bedGraph/$filename"
                else if( filename.indexOf("CpG_report" ) > 0 ) "stranded_CpG_report/$filename"
                else "methylation_calls/$filename"
            }

        input:
        set val(name), file(bam) from ch_bam_dedup_for_bismark_methXtract
        file index from ch_bismark_index_for_bismark_methXtract.collect()

        output:
        set val(name), file("*splitting_report.txt") into ch_bismark_splitting_report_for_bismark_report, ch_bismark_splitting_report_for_bismark_summary, ch_bismark_splitting_report_for_multiqc
        set val(name), file("*.M-bias.txt") into ch_bismark_mbias_for_bismark_report, ch_bismark_mbias_for_bismark_summary, ch_bismark_mbias_for_multiqc
        file '*.{png,gz}'

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        cytosine_report = params.cytosine_report ? "--cytosine_report --genome_folder ${index} " : ''
        meth_cutoff = params.meth_cutoff ? "--cutoff ${params.meth_cutoff}" : ''
        multicore = ''
        if( task.cpus ){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 3) as int
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }
        buffer = ''
        if( task.memory ){
            mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
            // only set if we have more than 6GB available
            if( mbuffer.compareTo(4.GB) == 1 ){
              buffer = "--buffer_size ${mbuffer.toGiga()}G"
            }
        }
        if(params.single_end) {
            """
            bismark_methylation_extractor $comprehensive $meth_cutoff \\
                $multicore $buffer $cytosine_report \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -s \\
                --report \\
                $bam
            """
        } else {
            """
            bismark_methylation_extractor $comprehensive $meth_cutoff \\
                $multicore $buffer $cytosine_report \\
                --ignore_r2 2 \\
                --ignore_3prime_r2 2 \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -p \\
                --no_overlap \\
                --report \\
                $bam
            """
        }
    }

    ch_bismark_align_log_for_bismark_report
     .join(ch_bismark_dedup_log_for_bismark_report)
     .join(ch_bismark_splitting_report_for_bismark_report)
     .join(ch_bismark_mbias_for_bismark_report)
     .set{ ch_bismark_logs_for_bismark_report }


    /*
     * STEP 6 - Bismark Sample Report
     */
    process bismark_report {
        tag "$name"
        publishDir "${params.outdir}/bismark_reports", mode: params.publish_dir_mode

        input:
        set val(name), file(align_log), file(dedup_log), file(splitting_report), file(mbias) from ch_bismark_logs_for_bismark_report

        output:
        file '*{html,txt}' into ch_bismark_reports_results_for_multiqc

        script:
        """
        bismark2report \\
            --alignment_report $align_log \\
            --dedup_report $dedup_log \\
            --splitting_report $splitting_report \\
            --mbias_report $mbias
        """
    }

    /*
     * STEP 7 - Bismark Summary Report
     */
    process bismark_summary {
        publishDir "${params.outdir}/bismark_summary", mode: params.publish_dir_mode

        input:
        file ('*') from ch_bam_for_bismark_summary.collect()
        file ('*') from ch_bismark_align_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_dedup_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_splitting_report_for_bismark_summary.collect()
        file ('*') from ch_bismark_mbias_for_bismark_summary.collect()

        output:
        file '*{html,txt}' into ch_bismark_summary_results_for_multiqc

        script:
        """
        bismark2summary
        """
    }
} // End of bismark processing block