#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include {help} from "${projectDir}/bin/help/illumina"
include {PARAM_LOG} from "${projectDir}/bin/log/illumina"
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {SAMTOOLS_FAIDX} from "${projectDir}/modules/samtools/samtools_faidx"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates_removedup"
include {SAMTOOLS_STATS} from "${projectDir}/modules/samtools/samtools_stats_mmrsvd"
include {SMOOVE_CALL} from "${projectDir}/modules/smoove/smoove_call_germline"
include {BCFTOOLS_REHEAD_SORT as REHEAD_SORT_LUMPY;
         BCFTOOLS_REHEAD_SORT as REHEAD_SORT_DELLY;
         BCFTOOLS_REHEAD_SORT as REHEAD_SORT_CNV;
         BCFTOOLS_REHEAD_SORT as REHEAD_SORT_MANTA} from "${projectDir}/modules/bcftools/bcftools_rehead_sort"
include {MANTA_CALL} from "${projectDir}/modules/illumina/manta_germline"
include {DELLY_CALL_GERMLINE} from "${projectDir}/modules/delly/delly_call_germline"
include {DELLY_CNV_GERMLINE} from "${projectDir}/modules/delly/delly_cnv_germline"

include {GATK_INDEXFEATUREFILE;
         GATK_INDEXFEATUREFILE as GATK_INDEXFEATUREFILE_SNV} from "${projectDir}/modules/gatk/gatk_indexfeaturefile"

include {GATK_HAPLOTYPECALLER_INTERVAL} from "${projectDir}/modules/gatk/gatk_haplotypecaller_interval"
include {MAKE_VCF_LIST} from "${projectDir}/modules/utility_modules/make_vcf_list"
include {GATK_MERGEVCF_LIST} from "${projectDir}/modules/gatk/gatk_mergevcf_list"
include {GATK_COMBINEGVCFS} from "${projectDir}/modules/gatk/gatk_combinegvcfs"
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration"
include {GATK_MERGEVCF} from "${projectDir}/modules/gatk/gatk_mergevcf"

include {VEP_GERMLINE as VEP_GERMLINE_GATK;
         VEP_GERMLINE as VEP_GERMLINE_CNV} from "${projectDir}/modules/ensembl/varianteffectpredictor_germline_mouse"

include {BCFTOOLS_VCF_TO_BCF} from "${projectDir}/modules/bcftools/bcftools_vcf_to_bcf"
include {DUPHOLD as DUPHOLD_DELLY;
         DUPHOLD as DUPHOLD_LUMPY;
         DUPHOLD as DUPHOLD_MANTA} from "${projectDir}/modules/duphold/duphold"
include {BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_DELLY;
         BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_LUMPY;
         BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_MANTA;} from "${projectDir}/modules/bcftools/bcftools_duphold_filter"

include {SV_MERGE} from "${projectDir}/modules/r/illumina_sv_merge"
include {PYTHON_BEDPE_TO_VCF} from "${projectDir}/modules/python/python_bedpe_to_vcf"
include {SURVIVOR_VCF_TO_TABLE} from "${projectDir}/modules/survivor/survivor_vcf_to_table"
include {SURVIVOR_SUMMARY} from "${projectDir}/modules/survivor/survivor_summary"
include {SURVIVOR_TO_BED} from "${projectDir}/modules/survivor/survivor_to_bed"
include {SURVIVOR_BED_INTERSECT} from "${projectDir}/modules/survivor/survivor_bed_intersect"
include {SURVIVOR_ANNOTATION} from "${projectDir}/modules/survivor/survivor_annotation"
include {SURVIVOR_INEXON} from "${projectDir}/modules/survivor/survivor_inexon"

// log parameter info
PARAM_LOG()

workflow ILLUMINA {

    if (params.help){
       help()
        exit 0
    }
 
    ch_fastq1 = params.fastq1 ? Channel.fromPath(params.fastq1) : null
    ch_fastq2 = params.fastq2 ? Channel.fromPath(params.fastq2) : null
    ch_sampleID = params.sampleID ? Channel.value(params.sampleID) : null
    ch_bam = params.bam ? Channel.fromPath(params.bam) : null

    // Prepare reads channel

    if (params.fastq1 && !params.fastq2 && !params.bam) {
        fq_reads = ch_sampleID.concat(ch_fastq1)
                            .collect()
                            .map { it -> tuple(it[0], it[1])}
    }

    else if (params.fastq1 && params.fastq2 && !params.bam) {
        fq_reads = ch_sampleID.concat(ch_fastq1, ch_fastq2)
                            .collect()
                            .map { it -> tuple(it[0], tuple(it[1], it[2]))}
    }

    else if (params.bam && !params.fastq1 && !params.fastq2) {
        fq_reads = null
        pre_bam = ch_sampleID.concat(ch_bam)
                             .collect()
                             .map { it -> tuple(it[0], it[1])}
    } else {
        exit 1, "Both FASTQ and BAM inputs were specified. Use either FASTQ or BAM as workflow input."
    }

    // Index reference fasta
    faidx_input = ['primary_strain', params.ref_fa]
    SAMTOOLS_FAIDX(faidx_input)
    fasta_index = SAMTOOLS_FAIDX.out.fai.map{it -> [params.ref_fa, it[1]]}

    // ** Optional mapping steps when input are FASTQ files
    if (params.fastq1) {
        // Filter and trim reads
        FASTP(fq_reads)
        
        // Get read groups ID from FASTQ file
        READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")

        // Map reads to reference
        bwa_mem_input = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
        BWA_MEM(bwa_mem_input)

        // Sort and compress to BAM
        SAMTOOLS_SORT(BWA_MEM.out.sam, '-O bam', 'bam')

        ch_bam_undup = SAMTOOLS_SORT.out.sorted_file
    }
    else {
        ch_bam_undup = pre_bam
    }

    // Remove optical duplicates from alignment
    PICARD_MARKDUPLICATES(ch_bam_undup)

    // Quantify insert sizes
    SAMTOOLS_STATS(PICARD_MARKDUPLICATES.out.bam_and_index)

    // Call SV with SMOOVE
    SMOOVE_CALL(PICARD_MARKDUPLICATES.out.bam_and_index)
    REHEAD_SORT_LUMPY(SMOOVE_CALL.out.lumpy_vcf, "lumpy", fasta_index)

    // * Manta

    // Call SV with Manta
    MANTA_CALL(PICARD_MARKDUPLICATES.out.bam_and_index, fasta_index)
    REHEAD_SORT_MANTA(MANTA_CALL.out.manta_sv, "manta", fasta_index)

    // * Delly

    // Call SV with Delly
    DELLY_CALL_GERMLINE(PICARD_MARKDUPLICATES.out.bam_and_index, fasta_index)
    REHEAD_SORT_DELLY(DELLY_CALL_GERMLINE.out.delly_bcf, "delly_sv", fasta_index)
   
    // Call CNV with Delly
    DELLY_CNV_GERMLINE(PICARD_MARKDUPLICATES.out.bam_and_index, fasta_index)
    REHEAD_SORT_CNV(DELLY_CNV_GERMLINE.out.delly_bcf, "delly_cnv", fasta_index)
    GATK_INDEXFEATUREFILE(REHEAD_SORT_CNV.out.vcf_sort)
    VEP_GERMLINE_CNV(REHEAD_SORT_CNV.out.vcf_sort.join(GATK_INDEXFEATUREFILE.out.idx)) // THE OUTPUTS FROM THESE SHOULD BE SAVED AND ARE CURRENTLY NOT

    // Duphold
    DUPHOLD_DELLY(PICARD_MARKDUPLICATES.out.bam_and_index.join(REHEAD_SORT_DELLY.out.vcf_sort), fasta_index, 'delly_sv') 
    DUPHOLD_MANTA(PICARD_MARKDUPLICATES.out.bam_and_index.join(REHEAD_SORT_MANTA.out.vcf_sort), fasta_index, 'manta')

    // Duphold Filter
    BCFTOOLS_DUPHOLD_FILTER_DELLY(DUPHOLD_DELLY.out.vcf, 'delly_sv')
    BCFTOOLS_DUPHOLD_FILTER_MANTA(DUPHOLD_MANTA.out.vcf, 'manta')
    BCFTOOLS_DUPHOLD_FILTER_LUMPY(REHEAD_SORT_LUMPY.out.vcf_sort, 'lumpy')

    // * Merge callers and annotate results

    // Join VCFs together by sampleID and run SURVIVOR merge

    survivor_input = BCFTOOLS_DUPHOLD_FILTER_DELLY.out.vcf.join(BCFTOOLS_DUPHOLD_FILTER_LUMPY.out.vcf).join(BCFTOOLS_DUPHOLD_FILTER_MANTA.out.vcf)
                     .map { it -> tuple(it[0], tuple(it[1], it[2], it[3]))}
    SV_MERGE(survivor_input)
    PYTHON_BEDPE_TO_VCF(SV_MERGE.out.bedpe)

    SURVIVOR_VCF_TO_TABLE(PYTHON_BEDPE_TO_VCF.out.vcf)
    SURVIVOR_SUMMARY(PYTHON_BEDPE_TO_VCF.out.vcf)
    bed_prep_input = SURVIVOR_VCF_TO_TABLE.out.annotation.join(SURVIVOR_SUMMARY.out.csv)
    SURVIVOR_TO_BED(bed_prep_input)
    SURVIVOR_BED_INTERSECT(SURVIVOR_TO_BED.out.sv_beds)
    surv_annot_input = SURVIVOR_TO_BED.out.sv_beds.join(SURVIVOR_BED_INTERSECT.out.intersected_beds).join(SURVIVOR_SUMMARY.out.csv).join(SURVIVOR_VCF_TO_TABLE.out.annotation)
    SURVIVOR_ANNOTATION(surv_annot_input)
    surv_inexon_input = PYTHON_BEDPE_TO_VCF.out.vcf.join(SURVIVOR_BED_INTERSECT.out.intersected_exons)
    SURVIVOR_INEXON(surv_inexon_input)
}
