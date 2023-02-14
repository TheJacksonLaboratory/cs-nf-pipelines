#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/sv.nf"
include {param_log} from "${projectDir}/bin/log/sv.nf"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {SHORT_ALIGNMENT_MARKING} from "${projectDir}/modules/nygc-short-alignment-marking/short_alignment_marking"
include {PICARD_CLEANSAM} from "${projectDir}/modules/picard/picard_cleansam"
include {PICARD_FIX_MATE_INFORMATION} from "${projectDir}/modules/picard/picard_fix_mate_information"
include {PICARD_MARKDUPLICATES}	from "${projectDir}/modules/picard/picard_markduplicates"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"
include {CONPAIR_TUMOR_PILEUP} from "${projectDir}/modules/conpair/conpair_tumor_pileup"
include {CONPAIR_NORMAL_PILEUP} from "${projectDir}/modules/conpair/conpair_normal_pileup"
include {CONPAIR} from "${projectDir}/modules/conpair/conpair"
include {GATK_HAPLOTYPECALLER_SV_GERMLINE} from "${projectDir}/modules/gatk/gatk_haplotypecaller_sv_germline"
include {GATK_SORTVCF_GERMLINE as GATK_SORTVCF_GERMLINE;
         GATK_SORTVCF_GERMLINE as GATK_SORTVCF_GENOTYPE} from "${projectDir}/modules/gatk/gatk_sortvcf_germline"
include {GATK_GENOTYPE_GVCF} from "${projectDir}/modules/gatk/gatk_genotype_gvcf"
include {GATK_CNNSCORE_VARIANTS} from "${projectDir}/modules/gatk/gatk_cnnscorevariants"
include {GATK_FILTER_VARIANT_TRANCHES} from "${projectDir}/modules/gatk/gatk_filtervarianttranches"
include {GATK_VARIANTFILTRATION_AF} from "${projectDir}/modules/gatk/gatk_variantfiltration_af"
include {BCFTOOLS_GERMLINE_FILTER} from "${projectDir}/modules/bcftools/bcftools_germline_filter"
include {BCFTOOLS_SPLITMULTIALLELIC_REGIONS} from "${projectDir}/modules/bcftools/bcftools_split_multiallelic_regions"
include {VEP_GERMLINE} from "${projectDir}/modules/ensembl/varianteffectpredictor"
include {BCFTOOLS_REMOVESPANNING} from "${projectDir}/modules/bcftools/bcftools_remove_spanning"
include {COSMIC_ANNOTATION} from "${projectDir}/modules/cosmic/cosmic_annotation"
include {COSMIC_CANCER_RESISTANCE_MUTATION} from "${projectDir}/modules/cosmic/cosmic_add_cancer_resistance_mutations"
include {GERMLINE_VCF_FINALIZATION} from "${projectDir}/modules/utility_modules/germline_vcf_finalization"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {SNPSIFT_EXTRACT_AND_PARSE} from "${projectDir}/modules/utility_modules/parse_extracted_sv_table"
include {GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_NORMAL;
         GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_TUMOR} from "${projectDir}/modules/gatk/gatk_getsamplename"
include {GATK_MUTECT2} from "${projectDir}/modules/gatk/gatk_mutect2"
include {GATK_MERGEMUTECTSTATS} from "${projectDir}/modules/gatk/gatk_mergemutectstats"
include {GATK_FILTERMUECTCALLS} from "${projectDir}/modules/gatk/gatk_filtermutectcalls"
include {MANTA} from "${projectDir}/modules/illumina/manta"
include {STRELKA2} from "${projectDir}/modules/illumina/strelka2"
include {LANCET} from "${projectDir}/modules/nygenome/lancet"
include {GATK_SORTVCF as GATK_SORTVCF_MUTECT;
         GATK_SORTVCF as GATK_SORTVCF_LANCET;
         GATK_SORTVCF as GATK_SORTVCF_TOOLS;
         GATK_SORTVCF as GATK_SORTVCF_TOOLS_LANCET} from "${projectDir}/modules/gatk/gatk_sortvcf_somatic_tools"
include {GRIDSS_PREPROCESS} from "${projectDir}/modules/gridss/gridss_preprocess"
include {GRIDSS_ASSEMBLE} from "${projectDir}/modules/gridss/gridss_assemble"
include {GRIDSS_CALLING} from "${projectDir}/modules/gridss/gridss_calling"
include {GRIDSS_CHROM_FILTER} from "${projectDir}/modules/gridss/gridss_chrom_filter"
include {GRIDSS_SOMATIC_FILTER} from "${projectDir}/modules/gridss/gridss_somatic_filter"
include {SAMTOOLS_STATS_INSERTSIZE as SAMTOOLS_STATS_INSERTSIZE_NORMAL;
         SAMTOOLS_STATS_INSERTSIZE as SAMTOOLS_STATS_INSERTSIZE_TUMOR} from "${projectDir}/modules/samtools/samtools_stats_insertsize"
include {SAMTOOLS_FILTER_UNIQUE as SAMTOOLS_FILTER_UNIQUE_NORMAL;
         SAMTOOLS_FILTER_UNIQUE as SAMTOOLS_FILTER_UNIQUE_TUMOR} from "${projectDir}/modules/samtools/samtools_filter_unique_reads"
include {BICSEQ2_NORMALIZE as BICSEQ2_NORMALIZE_NORMAL;
         BICSEQ2_NORMALIZE as BICSEQ2_NORMALIZE_TUMOR} from "${projectDir}/modules/biqseq2/bicseq2_normalize"
include {BICSEQ2_SEG} from "${projectDir}/modules/biqseq2/bicseq2_seg"
include {SVABA} from "${projectDir}/modules/svaba/svaba"
include {LUMPY_SV} from "${projectDir}/modules/lumpy_sv/lumpy_sv"
include {MSISENSOR2_MSI} from "${projectDir}/modules/msisensor2/msisensor2"

include {RENAME_METADATA;
         RENAME_METADATA as RENAME_METADATA_LANCET} from "${projectDir}/modules/python/python_rename_metadata"
include {MERGE_PREP;
         MERGE_PREP as MERGE_PREP_LANCET} from "${projectDir}/modules/python/python_merge_prep"
include {RENAME_VCF;
         RENAME_VCF as RENAME_VCF_LANCET;} from "${projectDir}/modules/python/python_rename_vcf"
include {COMPRESS_INDEX_VCF;
         COMPRESS_INDEX_VCF as COMPRESS_INDEX_VCF_LANCET;
         COMPRESS_INDEX_VCF as COMPRESS_INDEX_VCF_REGION_LANCET} from "${projectDir}/modules/tabix/compress_vcf"
include {BCFTOOLS_SPLITMULTIALLELIC;
         BCFTOOLS_SPLITMULTIALLELIC as BCFTOOLS_SPLITMULTIALLELIC_LANCET} from "${projectDir}/modules/bcftools/bcftools_split_multiallelic"
include {SPLIT_MNV;
         SPLIT_MNV as SPLIT_MNV_LANCET} from "${projectDir}/modules/python/python_split_mnv"
include {REMOVE_CONTIG} from "${projectDir}/modules/python/python_remove_contig"

include {BCFTOOLS_MERGECALLERS;
         BCFTOOLS_MERGECALLERS as BCFTOOLS_MERGECALLERS_FINAL} from "${projectDir}/modules/bcftools/bcftools_merge_callers"
include {BEDTOOLS_STARTCANDIDATES} from "${projectDir}/modules/bedtools/bedtools_start_candidates"
include {GET_CANDIDATES} from "${projectDir}/modules/python/python_get_candidates"
include {VCF_TO_BED} from "${projectDir}/modules/python/python_vcf_to_bed"
include {LANCET_CONFIRM} from "${projectDir}/modules/nygenome/lancet_confirm"
include {COMPRESS_INDEX_VCF_REGION;
         COMPRESS_INDEX_VCF_REGION as COMPRESS_INDEX_VCF_ALL_CALLERS;
         COMPRESS_INDEX_VCF_REGION as COMPRESS_INDEX_VCF_MERGED} from "${projectDir}/modules/tabix/compress_vcf_region"
include {BCFTOOLS_INTERSECTVCFS} from "${projectDir}/modules/bcftools/bcftools_intersect_lancet_candidates"

include {MERGE_COLUMNS} from "${projectDir}/modules/python/python_merge_columns"
include {ADD_NYGC_ALLELE_COUNTS} from "${projectDir}/modules/python/python_add_nygc_allele_counts"
include {ADD_FINAL_ALLELE_COUNTS} from "${projectDir}/modules/python/python_add_final_allele_counts"
include {FILTER_PON} from "${projectDir}/modules/python/python_filter_pon"
include {FILTER_VCF} from "${projectDir}/modules/python/python_filter_vcf"
include {SNV_TO_MNV_FINAL_FILTER} from "${projectDir}/modules/python/python_snv_to_mnv_final_filter"

include {GATK_SORTVCF_SOMATIC} from "${projectDir}/modules/gatk/gatk_sortvcf_somatic_merge"
include {REORDER_VCF_COLUMNS} from "${projectDir}/modules/python/python_reorder_vcf_columns"

// help if needed
if (params.help){
    help()
    exit 0
}

// log paramiter info
param_log()

// main workflow
workflow SV {

    if (params.csv_input) {
        ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))

        ch_input_sample.map{it -> [it[0], it[2]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 1: Qual_Stat
    QUALITY_STATISTICS(read_ch)

    // Step 2: Get Read Group Information
    READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "gatk")

    // Step 3: BWA-MEM Alignment
    bwa_mem_mapping = QUALITY_STATISTICS.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
    BWA_MEM(bwa_mem_mapping)
    
    // Step 4: Sort mapped reads
    PICARD_SORTSAM(BWA_MEM.out.sam)

    // Step 5: Remove short mapping 'artifacts': https://github.com/nygenome/nygc-short-alignment-marking
    SHORT_ALIGNMENT_MARKING(PICARD_SORTSAM.out.bam)

    // Step 6: Clean BAM to set MAPQ = 0 when read is unmapped (issue introduced in step 5)
    PICARD_CLEANSAM(PICARD_SORTSAM.out.bam)

    // Step 7: Fix mate information (fix pair flags due to mapping adjustment in step 5)
    PICARD_FIX_MATE_INFORMATION(PICARD_CLEANSAM.out.cleaned_bam)

    // Step 8: Markduplicates
    PICARD_MARKDUPLICATES(PICARD_FIX_MATE_INFORMATION.out.fixed_mate_bam)

    // Step 9: Calculate BQSR
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)

    // Step 10: Apply BQSR
    apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)
    GATK_APPLYBQSR(apply_bqsr)

    // Step 12: Nextflow channel processing
    // https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L854

    GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai).join(meta_ch).branch{
        normal: it[3].status == 0
        tumor:  it[3].status == 1
    }.set{chr_bam_status}
    // re-join the sampleID to metadata information. Split normal and tumor samples into 2 different paths. 
    // Process tumor and normal BAMs seperately for conpair. For calling, use mapped/joined data. 

    // Adjust channels to all normal, all tumor organized by patient IDs. 
    ch_bam_normal_to_cross = chr_bam_status.normal.map{ id, bam, bai, meta -> [meta.patient, meta, bam, bai] }
    ch_bam_tumor_to_cross = chr_bam_status.tumor.map{ id, bam, bai, meta -> [meta.patient, meta, bam, bai] }

    GATK_GETSAMPLENAME_NORMAL(ch_bam_normal_to_cross)
    GATK_GETSAMPLENAME_TUMOR(ch_bam_tumor_to_cross)

    ch_bam_normal_to_cross = ch_bam_normal_to_cross.join(GATK_GETSAMPLENAME_NORMAL.out.sample_name)
    ch_bam_tumor_to_cross = ch_bam_tumor_to_cross.join(GATK_GETSAMPLENAME_TUMOR.out.sample_name)

    // Cross all normal and tumor by patient ID. 
    ch_cram_variant_calling_pair = ch_bam_normal_to_cross.cross(ch_bam_tumor_to_cross)
        .map { normal, tumor ->
            def meta = [:]
            meta.patient    = normal[0]
            meta.normal_id  = normal[1].sample
            meta.tumor_id   = tumor[1].sample
            meta.sex        = normal[1].sex
            meta.id         = "${meta.tumor_id}_vs_${meta.normal_id}".toString()

            [normal[0], meta, normal[2], normal[3], normal[4], tumor[2], tumor[3], tumor[4]]
        }
        // normal[0] is patient ID, 
        // normal[1] and tumor[1] are meta info
        // normal[2] is normal bam, normal[3] is bai, normal[4] is read group ID. 
        // tumor[2] is bam, tumor[3] is bai, tumor[4] is read group ID.

    // Step 13: Conpair pileup for T/N
    CONPAIR_NORMAL_PILEUP(chr_bam_status.normal)
    CONPAIR_TUMOR_PILEUP(chr_bam_status.tumor)

    // output channel manipulation and cross/join
    conpair_normal_to_cross = CONPAIR_NORMAL_PILEUP.out.normal_pileup.map{ id, pileup, meta -> [meta.patient, meta, pileup] }
    conpair_tumor_to_cross = CONPAIR_TUMOR_PILEUP.out.tumor_pileup.map{ id, pileup, meta -> [meta.patient, meta, pileup] }

    conpair_input = conpair_normal_to_cross.cross(conpair_tumor_to_cross)
        .map { normal, tumor ->
            def meta = [:]
            meta.patient    = normal[0]
            meta.normal_id  = normal[1].sample
            meta.tumor_id   = tumor[1].sample
            meta.sex        = normal[1].sex
            meta.id         = "${meta.tumor_id}_vs_${meta.normal_id}".toString()

            [meta, normal[2], tumor[2]]
        }
        // normal[2] is normal pileup, tumor[2] is tumor pileup. 

    // the above channel manipulations will require a test in multiple sample mode. 

    // Step 12: Conpair for T/N concordance: https://github.com/nygenome/conpair
    CONPAIR(conpair_input)
    // NOTE: NEED HIGH COVERAGE TO TEST. 

    // Step 13: Germline Calling and annotation
    
    // Find the paths of all `scattered.interval_list` files, and make tuples with an index value. 
    // This is used for for HaplotypeCaller variant regions and GenotypeGVCF
    
    // Loads paths from a directory, collects into a list, sorts, 
    // maps to a set with indicies, flattens the map [file, index, file, index ...], 
    // collates the flattened map into pairs, then remaps to the pairs to tuple

    intervals = Channel.fromPath( params.chrom_intervals+'/*/scattered.interval_list' )
                .collect()
                .sort()
                .map { items -> items.withIndex() }
                .flatten()
                .collate(2)
                .map { item, idx -> tuple( item, idx + 1 ) }
    // https://stackoverflow.com/a/67084467/18557826

    // Applies scatter intervals from above to the BAM file channel prior to variant calling. 
    chrom_channel = ch_bam_normal_to_cross.combine(intervals)

    // Read a list of chromosome names from a parameter. These are provided to several tools. 
    chroms = Channel
        .fromPath("${params.chrom_contigs}")
        .splitText()
        .map{it -> it.trim()}

    // Get a list of primary chromosomes and exclude chrM (dropRight(1))
    chrom_list = chroms.collect().dropRight(1)
    chrom_list_noY = chrom_list.dropRight(1)
    
    // Variant calling. 
    GATK_HAPLOTYPECALLER_SV_GERMLINE(chrom_channel)

    // Applies gather to scattered haplotype calls.
    GATK_SORTVCF_GERMLINE(GATK_HAPLOTYPECALLER_SV_GERMLINE.out.vcf.groupTuple(), 'gvcf')

    // Applies scatter intervals from above to the merged file, and genotype.
    genotype_channel = GATK_SORTVCF_GERMLINE.out.vcf_idx.groupTuple().combine(intervals)
    // this will require a test in multiple sample mode. 

    GATK_GENOTYPE_GVCF(genotype_channel)
    GATK_CNNSCORE_VARIANTS(GATK_GENOTYPE_GVCF.out.vcf_idx)

    // Applies gather to genotyped/cnn called vcfs prior to tranche filtering. 
    GATK_SORTVCF_GENOTYPE(GATK_CNNSCORE_VARIANTS.out.vcf.groupTuple(), 'vcf')

    // Variant tranche filtering. 
    GATK_FILTER_VARIANT_TRANCHES(GATK_SORTVCF_GENOTYPE.out.vcf_idx)
    
    // Allele frequency and call refinement filtering. 
    GATK_VARIANTFILTRATION_AF(GATK_FILTER_VARIANT_TRANCHES.out.vcf_idx)
    BCFTOOLS_GERMLINE_FILTER(GATK_VARIANTFILTRATION_AF.out.vcf)

    // Germline annotation - Filtered
    // 1. SplitMultiAllelicRegions & compress & index
    BCFTOOLS_SPLITMULTIALLELIC_REGIONS(BCFTOOLS_GERMLINE_FILTER.out.vcf_idx, chrom_list_noY)
    // 2. vepPublicSvnIndel
    VEP_GERMLINE(BCFTOOLS_SPLITMULTIALLELIC_REGIONS.out.vcf_idx)
    // 3. RemoveSpanning
    BCFTOOLS_REMOVESPANNING(VEP_GERMLINE.out.vcf)
    // 4. AddCosmic
    COSMIC_ANNOTATION(BCFTOOLS_REMOVESPANNING.out.vcf)
    // 5. AddCancerResistanceMutations
    COSMIC_CANCER_RESISTANCE_MUTATION(COSMIC_ANNOTATION.out.vcf)
    // 6. AnnotateId & RenameCsqVcf
    GERMLINE_VCF_FINALIZATION(COSMIC_CANCER_RESISTANCE_MUTATION.out.vcf, 'filtered')

    SNPSIFT_EXTRACTFIELDS(GERMLINE_VCF_FINALIZATION.out.vcf)
    SNPSIFT_EXTRACT_AND_PARSE(SNPSIFT_EXTRACTFIELDS.out.temp)
    
    // NOTE: Annotation can be done on the GATK_VARIANTFILTRATION_AF.out.vcf_idx file
    //       The steps would need to be split with 'as' statements in the 'include' step, and then added here.

    // Step 14: Somatic Calling

    // Applies scatter intervals from above to the BQSR bam file
    somatic_calling_channel = ch_cram_variant_calling_pair.combine(intervals)

    // // Applies scatter intervals from above to the BQSR bam file
    // somatic_calling_channel = ch_cram_variant_calling_pair.combine(chroms)
    // // NOTE: The above code line will split by Mutect2 calling by indivdiaul chromosomes 'chroms'. 
    // //       Entire chromosomes are scattered. For WGS, this is computationally intensive. 
    // //       We changed to calling to be done based on the same intervals passed to the germline caller. 
    // //       These intervals are based on the 'noN' file make by BROAD/GATK. 
    // //       If complete chromosomes are requried, the above line of code can be uncommented. 

    // Mutect2
    // STEPS: Call on each chromosome / interval. 
    //        Prior to 'filtermutectcalls' vcfs must be merged (GATK best practice). 
    //             NOTE: The group and map statement ensures that VCFs are organzied by sampleID, and carry  and toolID is maintained through the process. 
    //        Prior to 'filtermutectcalls' "stats" files from mutect2 must be merged (GATK best practice).
    //        Merge vcfs and stats must be Nextflow joined prior to 'filtermutectcalls' to avoid samples being confounded. 

    GATK_MUTECT2(somatic_calling_channel)

    sort_merge_input_mutect2VCF = GATK_MUTECT2.out.vcf
                                  .groupTuple()
                                  .map {  sampleID, vcf, meta, normal, tumor, tool -> tuple( sampleID, vcf, meta.unique()[0], normal.unique()[0], tumor.unique()[0], tool.unique()[0] )  }
    
    GATK_SORTVCF_MUTECT(sort_merge_input_mutect2VCF)
    GATK_MERGEMUTECTSTATS(GATK_MUTECT2.out.stats.groupTuple())
    
    filter_mutect_input = GATK_SORTVCF_MUTECT.out.vcf_tbi.join(GATK_MERGEMUTECTSTATS.out.stats)

    GATK_FILTERMUECTCALLS(filter_mutect_input)
    // additional NYGC steps not used: add commands to VCF

    // Manta
    MANTA(ch_cram_variant_calling_pair)
    // additional NYGC steps not used: add commands to VCF
    // FilterNonpass is used in NYGC with `SelectVariants` and `--exclude-filtered`. Do we want hard filtering? 

    // Strelka2
    strekla2_input = ch_cram_variant_calling_pair.join(MANTA.out.manta_smallindel_vcf_tbi)
    STRELKA2(strekla2_input)
    // additional NYGC steps not used: add commands to VCF

    // Lancet
    // Generate a list of chromosome beds. This is generated in the same manner as the calling `intervals` variable above. 
    lancet_beds = Channel.fromPath( params.lancet_beds_directory+'/*.bed' )
                    .collect()
                    .sort()
                    .map { items -> items.withIndex() }
                    .flatten()
                    .collate(2)
                    .map { item, idx -> tuple( item, idx + 1 ) }
    // https://stackoverflow.com/a/67084467/18557826

    // Applies scatter intervals from above to the BQSR bam file
    lancet_calling_channel = ch_cram_variant_calling_pair.combine(lancet_beds)
    LANCET(lancet_calling_channel)
    // additional NYGC steps not used: add commands to VCF

    sort_merge_input_lancetVCF = LANCET.out.vcf
                                .groupTuple()
                                .map {  sampleID, vcf, meta, normal, tumor, tool -> tuple( sampleID, vcf, meta.unique()[0], normal.unique()[0], tumor.unique()[0], tool.unique()[0] )  }
    
    GATK_SORTVCF_LANCET(sort_merge_input_lancetVCF)

    // Gridss
    GRIDSS_PREPROCESS(ch_cram_variant_calling_pair)
    gridss_assemble_input = ch_cram_variant_calling_pair.join(GRIDSS_PREPROCESS.out.gridss_preproc)
    GRIDSS_ASSEMBLE(gridss_assemble_input)
    gridss_call_input = ch_cram_variant_calling_pair.join(GRIDSS_ASSEMBLE.out.gridss_assembly)
    GRIDSS_CALLING(gridss_call_input)
    GRIDSS_CHROM_FILTER(GRIDSS_CALLING.out.gridss_vcf, chrom_list)
    GRIDSS_SOMATIC_FILTER(GRIDSS_CHROM_FILTER.out.gridss_chrom_vcf, params.gridss_pon)
    // gridss somatic filter will need a higher coverage dataset for testing. 
    // additional NYGC steps not used: add commands to VCF 

    // BicSeq2
    SAMTOOLS_STATS_INSERTSIZE_NORMAL(ch_bam_normal_to_cross)
    SAMTOOLS_STATS_INSERTSIZE_TUMOR(ch_bam_tumor_to_cross)

    SAMTOOLS_FILTER_UNIQUE_NORMAL(ch_bam_normal_to_cross, chrom_list)
    SAMTOOLS_FILTER_UNIQUE_TUMOR(ch_bam_tumor_to_cross, chrom_list)

    biqseq_norm_input_normal = SAMTOOLS_FILTER_UNIQUE_NORMAL.out.uniq_seq.join(SAMTOOLS_STATS_INSERTSIZE_NORMAL.out.read_length_insert_size)
    // sampleID, individual_chr_seq_files, meta, read_ID, read_length, insert_size. 
    biqseq_norm_input_tumor = SAMTOOLS_FILTER_UNIQUE_TUMOR.out.uniq_seq.join(SAMTOOLS_STATS_INSERTSIZE_TUMOR.out.read_length_insert_size)
    // sampleID, individual_chr_seq_files, meta, read_ID, read_length, insert_size. 

    fasta_files = Channel.fromPath( file(params.ref_fa).parent + '/*_chr*' )
            .collect()
    // collect individual chr fasta files. These are located in the same directory as the main reference. 
    // if the extension of `name_chr#.fa` changes this match will break. 
    // Can this be made more flexible without requiring another input parameter?

    BICSEQ2_NORMALIZE_NORMAL(biqseq_norm_input_normal, fasta_files)
    BICSEQ2_NORMALIZE_TUMOR(biqseq_norm_input_tumor, fasta_files)
    // bicseq2 normalize will need a higher coverage dataset for testing. 

    bicseq2_seg_input = BICSEQ2_NORMALIZE_NORMAL.out.normalized_output
                            .join(BICSEQ2_NORMALIZE_TUMOR.out.normalized_output)
                            // .groupTuple()
                            .map{sampleID, norm_files, meta_1, norm_readID, tumor_files, meta_2, tumor_readID -> tuple(sampleID, norm_files, tumor_files, meta_1, norm_readID, tumor_readID)}
    // sampleID, individual_normal_norm_bin_files, individual_tumor_norm_bin_files, metadata, norm_readID, tumor_readID. 

    BICSEQ2_SEG(bicseq2_seg_input)
    
    // Svaba
    SVABA(ch_cram_variant_calling_pair)
    // NOTE: SVABA is not in calling_wkf.wdl or used in the 'merge' steps. If included an additonal step: RemoveContig must be run on SVABA vcfs. 

    // Lumpy
    LUMPY_SV(ch_cram_variant_calling_pair)
    // NOTE: LUMPY is not in calling_wkf.wdl or used in the 'merge' steps. 

    // Step 15: MSI
    MSISENSOR2_MSI(ch_bam_tumor_to_cross)
    
    /*
    The follow are the harmonized output channels for each tool: 

    Manta
    MANTA.out.manta_somaticsv_tbi

    Strelka_SV
    STRELKA2.out.strelka_snv_vcf_tbi

    Strelka_INDEL
    STRELKA2.out.strelka_indel_vcf_tbi

    Mutect2
    GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi

    Lancet
    GATK_SORTVCF_LANCET.out.lancet_vcf

    Gridss
    GRIDSS_SOMATIC_FILTER.out.gridss_filtered_bgz

    Bicseq2
    BICSEQ2_SEG.out.bicseq2_sv_calls
    */



    /*
        NOTE: 
        The next section of this workflow becomes highly complex. Files from each caller are passed through  
        a set of 'merge prep' steps. These steps apply various functions to manipulate the VCF header, 
        and also calls within the VCFs. Once the VCFs are prepared, a merge occurs. 
        Following the merge, non-exonic regions are parsed out, and calls in those regions are passed
        to Lancet for confirmation. Following this, confirmed calls are used as 'support' and merged back
        to the full caller call set. Additional manipulations are done on the VCF, and then the 'final' VCF
        is passed through to the annotation steps. Additional and different annotations are done on SV and CNV 
        calls. The steps are commented to faciliate understanding of what is being done. 
    */

    somatic_caller_concat = MANTA.out.manta_somaticsv_tbi.concat( STRELKA2.out.strelka_snv_vcf_tbi, 
                                                                  STRELKA2.out.strelka_indel_vcf_tbi, 
                                                                  GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 
                                                                  GATK_SORTVCF_LANCET.out.vcf_tbi )

    // Merge prep: 
    // 1. Rename VCF header to include tool name: 
    RENAME_METADATA(somatic_caller_concat)

    // 2. Order samples in VCF to 'normal', 'tumor' and prep for merge. 
    //    See script for list of changes applied to the VCF:
    MERGE_PREP(RENAME_METADATA.out.rename_metadata_vcf)

    // 3. Rename VCF header to specfied 'normal' and 'tumor' names, add tool prefix to sampleIDs. 
    RENAME_VCF(MERGE_PREP.out.merge_prep_vcf)

    // 4. Compress and Index VCF:
    COMPRESS_INDEX_VCF(RENAME_VCF.out.rename_vcf)

    // 5. Split out multi-allelic calls:
    BCFTOOLS_SPLITMULTIALLELIC(COMPRESS_INDEX_VCF.out.compressed_vcf_tbi)

    // 6. Split MNV calls:
    SPLIT_MNV(BCFTOOLS_SPLITMULTIALLELIC.out.vcf)

    // 7. Sort VCF:
    GATK_SORTVCF_TOOLS(SPLIT_MNV.out.split_mnv_vcf)

    callers_for_merge = GATK_SORTVCF_TOOLS.out.vcf_tbi
                        .groupTuple()
                        .map{sampleID, vcf, idx, meta, normal_sample, tumor_sample, tool_list -> tuple( sampleID, vcf, idx, meta.unique()[0] )  }
                        .combine(chrom_list_noY.flatten())
    // The above collects all callers on sampleID, then maps to avoid duplication of data and to drop the tool list, which is not needed anymore. 
    // Note that this could be done using the very 'by' in the groupTuple statement. However, the map is still required to remove the tool list. 


    // Merge Callers, Extract non-exonic calls and try to confirm those with Lancet, 
        // then prep confirmed calls for merged back to full merge set: 
    
    // ** Make all caller merge set, and compress and index:
    BCFTOOLS_MERGECALLERS(callers_for_merge)
    COMPRESS_INDEX_VCF_ALL_CALLERS(BCFTOOLS_MERGECALLERS.out.vcf)
    
    // ** Extract non-exonic, and try to confirm with Lancet. 
    // 1. Intersect with '-v' against a list of exonic regions. This step subsets calls to non-exonic regions. 
    BEDTOOLS_STARTCANDIDATES(BCFTOOLS_MERGECALLERS.out.vcf)

    // 2. Get candidates from intersected, using rules outlined in get_candidates.py (script docs provided by original dev).
    //    Compress and index the resulting VCF. 
    GET_CANDIDATES(BEDTOOLS_STARTCANDIDATES.out.vcf)
    COMPRESS_INDEX_VCF_REGION(GET_CANDIDATES.out.vcf)

    // 3. VCF to BED
    VCF_TO_BED(GET_CANDIDATES.out.vcf)
    
    // 4. Confirm extracted calls with Lancet: 
    //    Compress and index the resulting VCF.
    lancet_confirm_input = VCF_TO_BED.out.bed                         
                           .combine(ch_cram_variant_calling_pair, by: 0)
                           .map{sampleID, bed, meta, chrom, meta2, normal_bam, normal_bai, normal_name, tumor_bam, tumor_bai, tumor_name -> tuple( sampleID, bed, meta, normal_bam, normal_bai, normal_name, tumor_bam, tumor_bai, tumor_name, chrom )  }
    // The above combines output by sampleID with BAM files. Then maps to avoid duplication of data, and set input tuples for the steps that follow.  
    // Note that "combine" here, combines each output stream from VCF_TO_BED with ch_cram_variant_calling_pair, keeping the scattered chrom seperate. 

    LANCET_CONFIRM(lancet_confirm_input)
    COMPRESS_INDEX_VCF_REGION_LANCET(LANCET_CONFIRM.out.vcf)

    // 5. Intersect Lancet Confirm with candidate extractions. 
    candidate_lancet_intersect_input = COMPRESS_INDEX_VCF_REGION.out.compressed_vcf_tbi
                                       .join(COMPRESS_INDEX_VCF_REGION_LANCET.out.compressed_vcf_tbi, by: [0,6])
                                       .map{sampleID, chrom, vcf, tbi, meta, empty_name, empty_name2, vcf2, tbi2, meta2, normal_name, tumor_name -> tuple( sampleID, vcf, tbi, vcf2, tbi2, meta, normal_name, tumor_name, chrom )}
    // The above joins candidate VCF with Lancet Confirm VCF by sampleID and chrom. Then maps to avoid duplication of data, and set input tuples for the steps that follow.  
    // Note: A. The 'by' statement here, joins on sampleID and chrom, which correspond to index values 0 and 6 in the output tuples. 
    //       B. 'empty_name' is used here because 'normal_name' and 'tumor_name' are not required/used in the candidate steps. 
    //       C. 'normal_name' and 'tumor_name' are needed to match input tuple expectations for teh steps that follow. 

    BCFTOOLS_INTERSECTVCFS(candidate_lancet_intersect_input)

    lancet_confirm_mergePrep_input = BCFTOOLS_INTERSECTVCFS.out.vcf.map{sampleID, vcf, index, meta, normal_name, tumor_name -> tuple(sampleID, vcf, index, meta, normal_name, tumor_name, 'lancet_support')}
    // The above remaps the output tuple from BCFTOOLS_INTERSECTVCF to include the tool name 'lancet', which is needed for the steps that follow. 
    // 'lancet_support' is used to trigger `--support` in the MERGE_PREP_LANCET statement. Logic is present in RENAME_VCF_LANCET to set the header to 'lancet' rather than 'lancet_support'

    // ** Prep calls for merge back to all caller merge set. 
    // 1. Rename VCF header to include tool name: 
    RENAME_METADATA_LANCET(lancet_confirm_mergePrep_input)

    // 2. Order samples in VCF to 'normal', 'tumor' and prep for merge. 
    //    See script for list of changes applied to the VCF:
    //    This step is done as `--support`
    MERGE_PREP_LANCET(RENAME_METADATA_LANCET.out.rename_metadata_vcf)

    // 3. Rename VCF header to specfied 'normal' and 'tumor' names, add tool prefix to sampleIDs.
    RENAME_VCF_LANCET(MERGE_PREP_LANCET.out.merge_prep_vcf)

    // 4. Compress and Index VCF:
    COMPRESS_INDEX_VCF_LANCET(RENAME_VCF_LANCET.out.rename_vcf)

    // 5. Split out multi-allelic calls:
    BCFTOOLS_SPLITMULTIALLELIC_LANCET(COMPRESS_INDEX_VCF_LANCET.out.compressed_vcf_tbi)

    // 6. Split MNV calls:
    SPLIT_MNV_LANCET(BCFTOOLS_SPLITMULTIALLELIC_LANCET.out.vcf)

    // 7. Remove contig descriptions:
    REMOVE_CONTIG(SPLIT_MNV_LANCET.out.split_mnv_vcf)

    // 8. Sort VCF. 
    GATK_SORTVCF_TOOLS_LANCET(REMOVE_CONTIG.out.remove_contig_vcf)

    // ** Merge lancet confirmed back to all merged callers. Compress and index merged calls.  
    allCalls_lancetConfirm_merge_input = COMPRESS_INDEX_VCF_ALL_CALLERS.out.compressed_vcf_tbi
                                         .join(GATK_SORTVCF_TOOLS_LANCET.out.vcf_tbi, by: [0,6])
                                         .map{sampleID, chrom, vcf, tbi, meta, empty_name, empty_name2, vcf2, tbi2, meta2, normal_name, tumor_name -> tuple( sampleID, [vcf, vcf2], [tbi, tbi2], meta, chrom )}
    // BCFTOOLS_MERGE Requires an input tuple as follows: [val(sampleID), file(vcf), file(idx), val(meta), val(chrom)]
    // Join the output streams on sampleID and chrom, and then map to the require tuple structure. Note that [vcf, vcf2] makes a list that is understoon by the module. 

    BCFTOOLS_MERGECALLERS_FINAL(allCalls_lancetConfirm_merge_input)
    COMPRESS_INDEX_VCF_MERGED(BCFTOOLS_MERGECALLERS_FINAL.out.vcf)
    
    // ** Manipulation of VCF into final file to be passed to annotation modules. 
    // 1. Merge Columns. 
    //    See script merge_columns.py for the three features used in merge (script docs provided by original dev).
    MERGE_COLUMNS(COMPRESS_INDEX_VCF_MERGED.out.compressed_vcf_tbi)
    //  NOTE: !! !! !! The merge and re-arrangment require sanity checks. 
    
    // 2. Add Allele Count to VCF.
    //    "Runs pileup on tumor and normal bam files to compute allele counts for bi-allelic SNV and Indel variants in VCF file and adds pileup format columns to the VCF file.""
    addAlleleCounts_confirm_input = MERGE_COLUMNS.out.mergeColumn_vcf                         
                           .combine(ch_cram_variant_calling_pair, by: 0)
                           .map{sampleID, vcf, meta, chrom, meta2, normal_bam, normal_bai, normal_name, tumor_bam, tumor_bai, tumor_name -> tuple( sampleID, vcf, meta, normal_bam, normal_bai, tumor_bam, tumor_bai, chrom )  }
    ADD_NYGC_ALLELE_COUNTS(addAlleleCounts_confirm_input)

    // 3. Add Final Allele Counts to VCF
    ADD_FINAL_ALLELE_COUNTS(ADD_NYGC_ALLELE_COUNTS.out.vcf)
    // NOTE: !! !! !! Sanity check is needed. Based on the addition of '_indel' and '_sv' and '_support' to strelka2 and lancet call names. 

    // 4. Filter VCF based on PON
    FILTER_PON(ADD_FINAL_ALLELE_COUNTS.out.vcf)

    // 5. Filter VCF based on gnomad and "ALL_GRCh38_sites"
    FILTER_VCF(FILTER_PON.out.vcf)

    // 6. "SnvstomnvsCountsbasedfilterAnnotatehighconf" 
    //    Parses file and converts adjacent SNVs to MNVs if they have they match the MNV_ID and called_by fields.
    SNV_TO_MNV_FINAL_FILTER(FILTER_VCF.out.vcf)


    // ** Collect and Merge Chroms. 

    chrom_merge_input = SNV_TO_MNV_FINAL_FILTER.out.vcf
                        .groupTuple()
                        .map{sampleID, vcf, meta, chrom -> tuple( sampleID, vcf, meta.unique()[0] )  }
                        // Collect scattered chroms, remap to tuple without chrom names. 

    GATK_SORTVCF_SOMATIC(chrom_merge_input)
    REORDER_VCF_COLUMNS(GATK_SORTVCF_SOMATIC.out.vcf_idx)
    // output tuple = val(sampleID), path("*_mergePrep.vcf"), val(meta). 
    // meta = [patient:test, normal_id:test, tumor_id:test2, sex:XX, id:test2_vs_test] 
    //         This named list can be accessed in the script section prior to """ via calls like: meta.patient

    // Step NN: Get alignment and WGS metrics
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(GATK_APPLYBQSR.out.bam)
    PICARD_COLLECTWGSMETRICS(GATK_APPLYBQSR.out.bam)

}


// nextflow /projects/omics_share/meta/benchmarking/ngs-ops-nf-pipelines/main.nf -profile sumner --workflow sv --gen_org human --pubdir /projects/omics_share/meta/benchmarking/testing/sv -w /projects/omics_share/meta/benchmarking/testing/work --csv_input /projects/omics_share/meta/benchmarking/ngs-ops-nf-pipelines/sv_input.csv -resume -stub

// Function to extract information (meta data + file(s)) from csv file(s)
// https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L1084
def extract_csv(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    // Additional check of sample sheet:
    // 1. Each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def patient_sample_lane_combinations_in_samplesheet = []
    def sample2patient = [:]

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!sample2patient.containsKey(row.sample.toString())) {
                sample2patient[row.sample.toString()] = row.patient.toString()
            } else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
                log.error('The sample "' + row.sample.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sample.toString()] + '" in the sample sheet.')
                System.exit(1)
            }
        }

    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0

    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all++
            if (!(row.patient && row.sample)){
                log.error "Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'."
                System.exit(1)
            }
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no sex specified, sex is not considered
        // sex is only mandatory for somatic CNV
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 0) sample_count_normal++
        else sample_count_tumor++

        // join meta to fastq
        if (row.fastq_2) {
            meta.id         = "${row.patient}-${row.sample}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)

            meta.data_type  = 'fastq'

            meta.size       = 1 // default number of splitted fastq

            return [meta.id, meta, [fastq_1, fastq_2]]

        } else {
            log.error "Missing or unknown field in csv file header. Please check your samplesheet"
            System.exit(1)
        }
    }
}