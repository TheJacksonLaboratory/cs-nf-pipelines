#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include {help} from "${projectDir}/bin/help/illumina"
include {PARAM_LOG} from "${projectDir}/bin/log/illumina"
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {BWA_INDEX} from "${projectDir}/modules/bwa/bwa_index"
include {SAMTOOLS_FAIDX} from "${projectDir}/modules/samtools/samtools_faidx"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {GATK_MARK_DUPLICATES} from "${projectDir}/modules/gatk/gatk_mark_duplicates"
include {SAMTOOLS_STATS} from "${projectDir}/modules/samtools/samtools_stats"
include {LUMPY_PREP} from "${projectDir}/modules/lumpy/lumpy_prep"
include {PICARD_SORTSAM as LUMPY_SORTSAM;
         PICARD_SORTSAM as LUMPY_SORTSAM_DISCORDANT;
         PICARD_SORTSAM as LUMPY_SORTSAM_SPLIT} from "${projectDir}/modules/picard/picard_sortsam"
include {LUMPY_EXTRACT_SPLITS} from "${projectDir}/modules/lumpy/lumpy_extract_splits"
include {LUMPY_CALL_SV} from "${projectDir}/modules/lumpy/lumpy_call_sv"
include {REHEADER_VCF as REHEADER_LUMPY;
         REHEADER_VCF as REHEADER_DELLY} from "${projectDir}/modules/utility_modules/reheader_vcf"
include {MANTA_CALL} from "${projectDir}/modules/manta/manta_call"
include {DELLY_CALL} from "${projectDir}/modules/delly/delly_call"
include {DELLY_POST_PROCESS} from "${projectDir}/modules/delly/delly_post_process"
include {SURVIVOR_MERGE} from "${projectDir}/modules/survivor/survivor_merge"
include {SURVIVOR_VCF_TO_TABLE} from "${projectDir}/modules/survivor/survivor_vcf_to_table"
include {SURVIVOR_SUMMARY} from "${projectDir}/modules/survivor/survivor_summary"
include {SURVIVOR_TO_BED} from "${projectDir}/modules/survivor/survivor_to_bed"
include {SURVIVOR_BED_INTERSECT} from "${projectDir}/modules/survivor/survivor_bed_intersect"
include {SURVIVOR_ANNOTATION} from "${projectDir}/modules/survivor/survivor_annotation"
include {SURVIVOR_INEXON} from "${projectDir}/modules/survivor/survivor_inexon"

workflow ILLUMINA {

    if (params.help){
       help()
        exit 0
    }

    ch_fasta = params.fasta ? Channel.fromPath(params.fasta) : null
    ch_bwa_index = params.bwa_index ? Channel.fromPath(params.bwa_index) : null
    ch_fastq1 = params.fastq1 ? Channel.fromPath(params.fastq1) : null
    ch_fastq2 = params.fastq2 ? Channel.fromPath(params.fastq2) : null
    ch_sampleID = params.sampleID ? Channel.value(params.sampleID) : null
    ch_bam = params.bam ? Channel.fromPath(params.bam) : null
    ch_bwa_index = params.bwa_index ? Channel.fromPath(params.bwa_index) : null

    // Prepare reads channel

    if (params.fastq1 && !params.fastq2 && !params.bam && !params.sample_folder) {
        fq_reads = ch_sampleID.concat(ch_fastq1)
                            .collect()
                            .map { it -> tuple(it[0], it[1])}
    }

    else if (params.fastq1 && params.fastq2 && !params.bam && !params.sample_folder) {
        fq_reads = ch_sampleID.concat(ch_fastq1, ch_fastq2)
                            .collect()
                            .map { it -> tuple(it[0], tuple(it[1], it[2]))}
    }

    else if (params.bam && !params.fastq1 && !params.fastq2 && !params.sample_folder) {
        fq_reads = null
        pre_bam = ch_sampleID.concat(ch_bam)
                             .collect()
                             .map { it -> tuple(it[0], it[1])}
    }

    PARAM_LOG()

    // Index reference fasta
    SAMTOOLS_FAIDX(ch_fasta)
    
    // ** Optional mapping steps when input are FASTQ files
    if (params.fastq1) {
        
        // Generate reference index if neccesary
        if(!params.bwa_index) {
            BWA_INDEX(ch_fasta)
            ch_bwa_index = BWA_INDEX.out.bwa_index
        }
        else {
            ch_bwa_index = file("${params.bwa_index}.*")
        }

        // Filter and trim reads
        FASTP(fq_reads)
        
        // Get read groups ID from FASTQ file
        READ_GROUPS(FASTP.out.trimmed_fastq)

        // Map reads to reference
        //bwa_mem_input = FILTER_TRIM.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
        bwa_mem_input = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
        BWA_MEM(bwa_mem_input, SAMTOOLS_FAIDX.out.fasta_fai, ch_bwa_index)

        // Sort and compress to BAM
        SAMTOOLS_SORT(BWA_MEM.out.sam)

        ch_bam_undup = SAMTOOLS_SORT.out.bam
    }
    else {
        ch_bam_undup = pre_bam
    }

    // Remove optical duplicates from alignment
    GATK_MARK_DUPLICATES(ch_bam_undup)

    // Quantify insert sizes
    SAMTOOLS_STATS(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Prep BAM for Lumpy (Map clipped reads, read group info, extract discordant alignments)
    LUMPY_PREP(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Sort prepped LUMPY bams
    LUMPY_SORTSAM(LUMPY_PREP.out.bam_bwa_lumpy, "alignBWA_lumpySort")
    LUMPY_SORTSAM_DISCORDANT(LUMPY_PREP.out.dis_unsorted_bam, "alignBWA_lumpySort_discordant")

    // Extract split reads
    LUMPY_EXTRACT_SPLITS(LUMPY_PREP.out.bam_bwa_lumpy)
    LUMPY_SORTSAM_SPLIT(LUMPY_EXTRACT_SPLITS.out.bam_bwa_lumpy, "alignBWA_lumpySort_splits")

    // Call SV with Lumpy
    lumpy_input = LUMPY_SORTSAM.out.sorted_bam.join(LUMPY_SORTSAM_SPLIT.out.sorted_bam).join(LUMPY_SORTSAM_DISCORDANT.out.sorted_bam)
    LUMPY_CALL_SV(lumpy_input)
    REHEADER_LUMPY(LUMPY_CALL_SV.out.lumpy_vcf, "lumpy")

    // * Manta

    // Call SV with Manta
    MANTA_CALL(GATK_MARK_DUPLICATES.out.bam_and_index, SAMTOOLS_FAIDX.out.fasta_fai)

    // * Delly

    // Call SV with Delly
    DELLY_CALL(GATK_MARK_DUPLICATES.out.bam_and_index, SAMTOOLS_FAIDX.out.fasta_fai)
    
    // Convert Delly BCF to VCF and reheader with sample name
    DELLY_POST_PROCESS(DELLY_CALL.out.delly_bcf)
    REHEADER_DELLY(DELLY_POST_PROCESS.out.delly_vcf, "delly")

    // * Merge callers and annotate results

    // Join VCFs together by sampleID and run SURVIVOR merge

    survivor_input = REHEADER_DELLY.out.vcf_rehead.join(REHEADER_LUMPY.out.vcf_rehead).join(MANTA_CALL.out.manta_sv)
                     .map { it -> tuple(it[0], tuple(it[1], it[2], it[3]))}
    SURVIVOR_MERGE(survivor_input)
    SURVIVOR_VCF_TO_TABLE(SURVIVOR_MERGE.out.vcf)
    SURVIVOR_SUMMARY(SURVIVOR_MERGE.out.vcf)

    bed_prep_input = SURVIVOR_VCF_TO_TABLE.out.annotation.join(SURVIVOR_SUMMARY.out.csv)
    SURVIVOR_TO_BED(bed_prep_input)
    SURVIVOR_BED_INTERSECT(SURVIVOR_TO_BED.out.sv_beds)
    surv_annot_input = SURVIVOR_TO_BED.out.sv_beds.join(SURVIVOR_BED_INTERSECT.out.intersected_beds).join(SURVIVOR_SUMMARY.out.csv).join(SURVIVOR_VCF_TO_TABLE.out.annotation)
    SURVIVOR_ANNOTATION(surv_annot_input)
    surv_inexon_input = SURVIVOR_MERGE.out.vcf.join(SURVIVOR_BED_INTERSECT.out.intersected_exons)
    SURVIVOR_INEXON(surv_inexon_input)
}