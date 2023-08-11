#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/pta.nf"
include {param_log} from "${projectDir}/bin/log/pta.nf"
include {CONCATENATE_PTA_FASTQ} from "${projectDir}/subworkflows/concatenate_pta_fastq"
include {JAX_TRIMMER} from "${projectDir}/modules/utility_modules/jax_trimmer"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {XENOME_CLASSIFY} from "${projectDir}/modules/xenome/xenome"
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

include {CONPAIR_PILEUP as CONPAIR_TUMOR_PILEUP;
         CONPAIR_PILEUP as CONPAIR_NORMAL_PILEUP} from "${projectDir}/modules/conpair/conpair_pileup"
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
include {VEP_GERMLINE} from "${projectDir}/modules/ensembl/varianteffectpredictor_germline"
include {BCFTOOLS_REMOVESPANNING} from "${projectDir}/modules/bcftools/bcftools_remove_spanning"
include {COSMIC_ANNOTATION} from "${projectDir}/modules/cosmic/cosmic_annotation"
include {COSMIC_CANCER_RESISTANCE_MUTATION_GERMLINE} from "${projectDir}/modules/cosmic/cosmic_add_cancer_resistance_mutations_germline"
include {GERMLINE_VCF_FINALIZATION} from "${projectDir}/modules/python/python_germline_vcf_finalization"
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
include {GRIPSS_SOMATIC_FILTER} from "${projectDir}/modules/gridss/gripss_somatic_filter"
include {SAMTOOLS_STATS_INSERTSIZE as SAMTOOLS_STATS_INSERTSIZE_NORMAL;
         SAMTOOLS_STATS_INSERTSIZE as SAMTOOLS_STATS_INSERTSIZE_TUMOR} from "${projectDir}/modules/samtools/samtools_stats_insertsize"
include {SAMTOOLS_FILTER_UNIQUE as SAMTOOLS_FILTER_UNIQUE_NORMAL;
         SAMTOOLS_FILTER_UNIQUE as SAMTOOLS_FILTER_UNIQUE_TUMOR} from "${projectDir}/modules/samtools/samtools_filter_unique_reads"
include {BICSEQ2_NORMALIZE as BICSEQ2_NORMALIZE_NORMAL;
         BICSEQ2_NORMALIZE as BICSEQ2_NORMALIZE_TUMOR} from "${projectDir}/modules/biqseq2/bicseq2_normalize"
include {BICSEQ2_SEG} from "${projectDir}/modules/biqseq2/bicseq2_seg"
include {BICSEQ2_SEG_UNPAIRED} from "${projectDir}/modules/biqseq2/bicseq2_seg_unpaired"
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
include {COMPRESS_INDEX_MERGED_VCF} from "${projectDir}/modules/tabix/compress_merged_vcf"
include {VEP_SOMATIC} from "${projectDir}/modules/ensembl/varianteffectpredictor_somatic"
include {COSMIC_ANNOTATION_SOMATIC} from "${projectDir}/modules/cosmic/cosmic_annotation_somatic"
include {COSMIC_CANCER_RESISTANCE_MUTATION_SOMATIC} from "${projectDir}/modules/cosmic/cosmic_add_cancer_resistance_mutations_somatic"
include {SOMATIC_VCF_FINALIZATION} from "${projectDir}/modules/python/python_somatic_vcf_finalization"
include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP_GERMLINE;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP_SOMATIC} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {ANNOTATE_BICSEQ2_CNV} from "${projectDir}/modules/r/annotate_bicseq2_cnv"
include {MERGE_SV} from "${projectDir}/modules/r/merge_sv"
include {ANNOTATE_SV;
         ANNOTATE_SV as ANNOTATE_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/annotate_sv"
include {ANNOTATE_GENES_SV;
         ANNOTATE_GENES_SV as ANNOTATE_GENES_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/annotate_genes_sv"
include {ANNOTATE_SV_WITH_CNV;
         ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL} from "${projectDir}/modules/r/annotate_sv_with_cnv"
include {FILTER_BEDPE;
         FILTER_BEDPE as FILTER_BEDPE_SUPPLEMENTAL} from "${projectDir}/modules/r/filter_bedpe"

include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"


// help if needed
if (params.help){
    help()
    exit 0
}

// log paramiter info
param_log()

// main workflow
workflow PTA {

    if (params.csv_input) {
        ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
        // Concat local Fastq files from CSV input if required.
            CONCATENATE_PTA_FASTQ(ch_input_sample)
            CONCATENATE_PTA_FASTQ.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
            CONCATENATE_PTA_FASTQ.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // ** Step 1: Qual_Stat
    JAX_TRIMMER(read_ch)
    
    FASTQC(JAX_TRIMMER.out.trimmed_fastq)

    // ** Step 2: Get Read Group Information
    READ_GROUPS(JAX_TRIMMER.out.trimmed_fastq, "gatk")

    // PDX CASES TO ADD AND VALIDATE: 
    // Normal samples should PASS the PDX step. 

    // ** Step 2a: Xenome if PDX data used.
    ch_XENOME_CLASSIFY_multiqc = Channel.empty() //optional log file. 
    if (params.pdx){
        // Xenome Classification
        XENOME_CLASSIFY(JAX_TRIMMER.out.trimmed_fastq)
        ch_XENOME_CLASSIFY_multiqc = XENOME_CLASSIFY.out.xenome_stats // set log file for multiqc

        bwa_mem_mapping = XENOME_CLASSIFY.out.xenome_human_fastq.join(READ_GROUPS.out.read_groups)

    } else { 
        bwa_mem_mapping = JAX_TRIMMER.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
    }

    // ** Step 3: BWA-MEM Alignment
    BWA_MEM(bwa_mem_mapping)
    
    // ** Step 4: Sort mapped reads
    PICARD_SORTSAM(BWA_MEM.out.sam)

    // ** Step 5: Remove short mapping 'artifacts': https://github.com/nygenome/nygc-short-alignment-marking
    SHORT_ALIGNMENT_MARKING(PICARD_SORTSAM.out.bam)

    // ** Step 6: Clean BAM to set MAPQ = 0 when read is unmapped (issue introduced in step 5)
    PICARD_CLEANSAM(PICARD_SORTSAM.out.bam)

    // ** Step 7: Fix mate information (fix pair flags due to mapping adjustment in step 5)
    PICARD_FIX_MATE_INFORMATION(PICARD_CLEANSAM.out.cleaned_bam)

    // ** Step 8: Markduplicates
    PICARD_MARKDUPLICATES(PICARD_FIX_MATE_INFORMATION.out.fixed_mate_bam)

    // ** Step 9: Calculate BQSR
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)

    // ** Step 10: Apply BQSR
    apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)
    GATK_APPLYBQSR(apply_bqsr)

    // Step 12: Nextflow channel processing
    // https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L854

    GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai).join(meta_ch).branch{
        normal: it[3].status == 0
        tumor:  it[3].status == 1
    }.set{ch_final_bam}
    // re-join the sampleID to metadata information. Split normal and tumor samples into 2 different paths. 
    // Process tumor and normal BAMs seperately for conpair. For calling, use mapped and crossed data. 

    // ** Get alignment and WGS metrics
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(GATK_APPLYBQSR.out.bam)
    PICARD_COLLECTWGSMETRICS(GATK_APPLYBQSR.out.bam)

    // ** NEXTFLOW OPERATORS::: Establish channels with sample pairs and individual input objects for downstream calling

    // get sample names, and join to bams. 
    GATK_GETSAMPLENAME_NORMAL(ch_final_bam.normal.map{ id, bam, bai, meta -> [id, meta, bam, bai] })
    GATK_GETSAMPLENAME_TUMOR(ch_final_bam.tumor.map{ id, bam, bai, meta -> [id, meta, bam, bai] })

    ch_normal_to_cross = ch_final_bam.normal.join(GATK_GETSAMPLENAME_NORMAL.out.sample_name).map{ id, bam, bai, meta, readID -> [meta.patient, meta, bam, bai, readID] }
    ch_tumor_to_cross  = ch_final_bam.tumor.join(GATK_GETSAMPLENAME_TUMOR.out.sample_name).map{ id, bam, bai, meta, readID -> [meta.patient, meta, bam, bai, readID] }
    
    /* 
    The above map statements adjusts channels for normal, tumor samples to organize them by patient IDs. 
    A common key ID is needed to cross tumor by normal samples when multiples of each patient are run together.

    NOTE!!!! that if a common patient key is then used as 'sampleID' across samples, 
    downstream results will have name collisions and results will be overwritten in muliple tumor to normal mappings. 

    e.g., patient: foo; tumor1 = bar, tumor2 = baz; normal = fizz. 
    
    Common key: sampleID == patient, results = foo.calls.vcf for both bar--fizz and baz--fizz, and results are mangled. 
    
    Unique key: sampleID == patient--tumor--normal, results == foo--bar--fizz.calls.vcf & foo--baz--fizz.calls.vcf. Results OK. 

    Therefore, the above ch_*_to_cross should ONLY be used in crossing samples. 
    A different channel is made below for cases when needed in callers. 
    */ 

    // Cross all normal and tumor by common patient ID. 
    ch_paired_samples = ch_normal_to_cross.cross(ch_tumor_to_cross)
        .map { normal, tumor ->
            def meta = [:]
            meta.patient    = normal[0]
            meta.normal_id  = normal[1].sampleID
            meta.tumor_id   = tumor[1].sampleID
            meta.sex        = normal[1].sex
            meta.id         = "${meta.patient}--${meta.tumor_id}--${meta.normal_id}".toString()

            [meta.id, meta, normal[2], normal[3], normal[4], tumor[2], tumor[3], tumor[4]]
        }
        /*
            normal[0] is patient ID, 
            normal[1] and tumor[1] are meta info
            normal[2] is normal bam, normal[3] is bai, normal[4] is read group ID. 
            tumor[2] is bam, tumor[3] is bai, tumor[4] is read group ID.
        */

    // Restore un-paired tumor samples, and add NA12878 as pairing in those cases
    ch_paired_samples = ch_tumor_to_cross
        .mix(ch_paired_samples)
        .map{it -> [it[1].patient, it[1], it[2], it[3], it[4]]}.groupTuple().filter{it[2].size() == 1} 
                    // it[0] = sampleID, it[1] = meta, it[2] = bam, it[3] = bai, it[4] = sampleReadID. 
                    // unknown group size, no 'size' statement can be used in groupTuple
        .map{tumor -> 
        def meta = [:]
            meta.patient    = tumor[1][0].patient
            meta.normal_id  = 'NA12878'
            meta.tumor_id   = tumor[1][0].sampleID
            meta.sex        = tumor[1][0].sex
            meta.id         = "${meta.patient}--${meta.tumor_id}--${meta.normal_id}".toString()
        
            [meta.id, meta, params.na12878_bam, params.na12878_bai, params.na12878_sampleName, tumor[2][0], tumor[3][0], tumor[4][0]]
        }
        .mix(ch_paired_samples)
        
    /* SAMPLE PAIRING CASES AND NOTES:
        1. Paired only for all samples: Managed by the intial cross statement. 
        2. Tumor only provide for all samples: Managed by the intial cross and subsequent remapping of non-crossed samples. 
        3. Some samples have a pair and others do not: Managed by the intial cross, and subsequent remapping of non-crossed samples.
        4. Multiple tumors per normal, or multiple normals per tumor, or a mixture of this: See note below. 

        Notes: 
        The cross statement manages one normal to many tumors, many normals to one tumor, and many normals to many tumors. 
        E.g.,:  
            [foo, [patient:foo, normal_id:n_baz, tumor_id:t_bar, sex:XX, id:t_bar_vs_n_baz], ....bam, ....bai, <RG1>, ....bam, ....bai, <RG2>]
            [foo, [patient:foo, normal_id:n_baz, tumor_id:t_qux, sex:XX, id:t_qux_vs_n_baz], ....bam, ....bai, <RG1>, ....bam, ....bai, <RG2>]
            ...
 
        When samples are provided without a pair, they will not be paired in the cross statment and dropped from the first: 'ch_paired_samples' instantiation. 
        To recover un-paired tumors, pair them with NA12878 and pass them with paired samples downstream,
        the group of all tumor samples: 'ch_tumor_to_cross' is mixed with the paired sample: 'ch_paired_samples'.
        Cases where tumor samples were paired are then filtered. Tumors that were paired in the cross will appear > 2 times in the mix results, and are removed via the 'filter it[2].size()==1' statement. 
        The resulting tumor-only samples are mapped into the format seen in the cross statement, with NA12878 being added via parameters as the 'normal' sample. 
    */


    ch_ind_samples = ch_paired_samples
        .filter{it[4] != params.na12878_sampleName}
        .multiMap{it -> 
                normal: ["${it[1].patient}--${it[1].normal_id}".toString(), it[1], it[2], it[3], it[4]]
                tumor:  ["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]
                }
        ch_normal_samples = ch_ind_samples.normal.unique{it[0]}
        ch_tumor_samples  = ch_ind_samples.tumor.unique{it[0]}

    ch_tumor_only = ch_paired_samples
        .filter{it[4] == params.na12878_sampleName}
        .map{it -> ["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]}
        .unique{it[0]}

    ch_msisensor2_input = ch_paired_samples
        .map{["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]}
        .unique{it[0]}

    /*
        The above establishes channels needed for germline calling, bicseq2 and MSIsensor2. 
        Those steps require BAM, index and readID. 
        Here sampleID is reset to the original ID from the CSV parser which is: 'patient--sample'
        Note that NA12878 is filtered for germline and can be filtered for bicseq2 / conpair,
        CNA and sample comparision analysis may not make sense for that pairing. 
        All tumor samples are passed to MSIsensor2 as it runs in tumor-only mode. 
    */



    // ** Step 13: Conpair pileup for T/N true pairs. 
    //    Step not run on tumor-only samples. As contamination analysis is not biologcally relavent. 

    conpair_input = ch_paired_samples
        .filter{it[4] != params.na12878_sampleName}
        .multiMap{it -> 
                normal: [it[1].patient, "${it[1].normal_id}".toString(), it[2], it[3]]
                tumor:  [it[1].patient, "${it[1].tumor_id}".toString(), it[5], it[6]]
                }
    /* 
        Remap the paired samples to required normal/tumor inputs for conpair, and filter NA12878 paired samples. 
        it[1] = metadata, it[2] = normal BAM, it[3] = normal BAI. 
        it[5] = tumor BAM, it[6] = tumor BAI.
        Patient ID is used here because samples must be re-crossed after the pileup to match all tumors and normals. 
    */ 

    CONPAIR_NORMAL_PILEUP(conpair_input.normal.unique{it[2]}, 'normal')
    CONPAIR_TUMOR_PILEUP(conpair_input.tumor.unique{it[2]}, 'tumor') 

    conpair_input = CONPAIR_NORMAL_PILEUP.out.pileup.cross(CONPAIR_TUMOR_PILEUP.out.pileup)
        .map { normal, tumor -> [normal[0], "${normal[0]}--${tumor[1]}--${normal[1]}".toString(), normal[2], tumor[2]]
        }
        // normal[0] is patientID or 'sampleID', normal[2] is normal pileup, tumor[2] is tumor pileup. 

    CONPAIR(conpair_input)


    // ** Step 14: Germline Calling and Annotation
    
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
    interval_count = files( params.chrom_intervals+'/*/scattered.interval_list' ).size()
    // interval count is used in groupTuple size statements. 

    // Applies scatter intervals from above to the BAM file channel prior to variant calling. 
    chrom_channel = ch_normal_samples.combine(intervals).filter{it[4] != params.na12878_sampleName}

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
    GATK_SORTVCF_GERMLINE(GATK_HAPLOTYPECALLER_SV_GERMLINE.out.vcf.groupTuple(size: interval_count), 'gvcf')

    // Applies scatter intervals from above to the merged file, and genotype.
    genotype_channel = GATK_SORTVCF_GERMLINE.out.vcf_idx.combine(intervals)

    GATK_GENOTYPE_GVCF(genotype_channel)
    GATK_CNNSCORE_VARIANTS(GATK_GENOTYPE_GVCF.out.vcf_idx)

    // Applies gather to genotyped/cnn called vcfs prior to tranche filtering. 
    GATK_SORTVCF_GENOTYPE(GATK_CNNSCORE_VARIANTS.out.vcf.groupTuple(size: interval_count), 'vcf')

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
    // 5. AddCancerResistanceMutations and dbsnpIDs
    COSMIC_CANCER_RESISTANCE_MUTATION_GERMLINE(COSMIC_ANNOTATION.out.vcf)
    SNPSIFT_ANNOTATE_DBSNP_GERMLINE(COSMIC_CANCER_RESISTANCE_MUTATION_GERMLINE.out.vcf, params.dbSNP, params.dbSNP_index, 'intermediate')
    // 6. AnnotateId & RenameCsqVcf
    GERMLINE_VCF_FINALIZATION(SNPSIFT_ANNOTATE_DBSNP_GERMLINE.out.vcf, 'filtered')
  


    // ** Step 15: Somatic Calling

    // Applies scatter intervals from above to the BQSR bam file
    somatic_calling_channel = ch_paired_samples.combine(intervals)

    /* Applies scatter intervals from above to the BQSR bam file
        somatic_calling_channel = ch_paired_samples.combine(chroms)
        NOTE: The above code line will split by Mutect2 calling by individual 'chroms'. 
            Entire chromosomes are scattered. For WGS, this is computationally intensive. 
            We changed to calling to be done based on the same intervals passed to the germline caller. 
            These intervals are based on the 'NoN' file made by BROAD/GATK. 
            If complete chromosomes are requried, the above line of code can be uncommented. 
       */

    // ** Mutect2 - SNP/InDEL Calling
    // STEPS: Call on each chromosome / interval. 
    //        Prior to 'filtermutectcalls' vcfs must be merged (GATK best practice). 
    //             NOTE: The group and map statement ensures that VCFs are organzied by sampleID, and carry  and toolID is maintained through the process. 
    //        Prior to 'filtermutectcalls' "stats" files from mutect2 must be merged (GATK best practice).
    //        Merge vcfs and stats must be Nextflow joined prior to 'filtermutectcalls' to avoid samples being confounded. 

    GATK_MUTECT2(somatic_calling_channel)

    sort_merge_input_mutect2VCF = GATK_MUTECT2.out.vcf
                                  .groupTuple(size: interval_count)
                                  .map {  sampleID, vcf, meta, normal, tumor, tool -> tuple( sampleID, vcf, meta.unique()[0], normal.unique()[0], tumor.unique()[0], tool.unique()[0] )  }
    
    GATK_SORTVCF_MUTECT(sort_merge_input_mutect2VCF)
    GATK_MERGEMUTECTSTATS(GATK_MUTECT2.out.stats.groupTuple(size: interval_count))
 
    filter_mutect_input = GATK_SORTVCF_MUTECT.out.vcf_tbi.join(GATK_MERGEMUTECTSTATS.out.stats)

    GATK_FILTERMUECTCALLS(filter_mutect_input)

    // ** Lancet - SNP/InDEL Calling
    // Generate a list of chromosome beds. This is generated in the same manner as the calling `intervals` variable above. 
    lancet_beds = Channel.fromPath( params.lancet_beds_directory+'/*.bed' )
                    .collect()
                    .sort()
                    .map { items -> items.withIndex() }
                    .flatten()
                    .collate(2)
                    .map { item, idx -> tuple( item, idx + 1 ) }
    // https://stackoverflow.com/a/67084467/18557826
    lancet_beds_count = files( params.lancet_beds_directory+'/*.bed' ).size()
    // bed file count is used in groupTuple size statements. 

    // Applies scatter intervals from above to the BQSR bam file
    lancet_calling_channel = ch_paired_samples.combine(lancet_beds)
    LANCET(lancet_calling_channel)

    sort_merge_input_lancetVCF = LANCET.out.vcf
                                .groupTuple(size: lancet_beds_count)
                                .map {  sampleID, vcf, meta, normal, tumor, tool -> tuple( sampleID, vcf, meta.unique()[0], normal.unique()[0], tumor.unique()[0], tool.unique()[0] )  }

    GATK_SORTVCF_LANCET(sort_merge_input_lancetVCF)

    // ** Manta - SV Calling
    MANTA(ch_paired_samples)
    // FilterNonpass can be used with `SelectVariants` and `--exclude-filtered`. However, hard filtering excluded for now. 

    // ** Strelka2 - SNP/InDEL Calling
    strekla2_input = ch_paired_samples.join(MANTA.out.manta_smallindel_vcf_tbi)
    STRELKA2(strekla2_input)

    // ** Gridss - SV Calling
    GRIDSS_PREPROCESS(ch_paired_samples)
    gridss_assemble_input = ch_paired_samples.join(GRIDSS_PREPROCESS.out.gridss_preproc)
    GRIDSS_ASSEMBLE(gridss_assemble_input)
    gridss_call_input = ch_paired_samples.join(GRIDSS_ASSEMBLE.out.gridss_assembly)
    GRIDSS_CALLING(gridss_call_input)
    GRIDSS_CHROM_FILTER(GRIDSS_CALLING.out.gridss_vcf, chrom_list)
    GRIPSS_SOMATIC_FILTER(GRIDSS_CHROM_FILTER.out.gridss_chrom_vcf)
    // NOTE: this filtering tool is hard coded for GRCh38 based on PON naming. 

    // ** BicSeq2 - CNV Calling
    /* This step does not run on unpaired samples. 
      CNV of tumor samples against an unrelated normal will produce spurious results. 
      NA12878 paired samples can be filtered from ch_normal_samples and ch_tumor_samples channels at their creation. 
    */
    SAMTOOLS_STATS_INSERTSIZE_NORMAL(ch_normal_samples)
    SAMTOOLS_STATS_INSERTSIZE_TUMOR(ch_tumor_samples.mix(ch_tumor_only))

    SAMTOOLS_FILTER_UNIQUE_NORMAL(ch_normal_samples, chrom_list)
    SAMTOOLS_FILTER_UNIQUE_TUMOR(ch_tumor_samples.mix(ch_tumor_only), chrom_list)

    biqseq_norm_input_normal = SAMTOOLS_FILTER_UNIQUE_NORMAL.out.uniq_seq.join(SAMTOOLS_STATS_INSERTSIZE_NORMAL.out.read_length_insert_size)
    // sampleID, individual_chr_seq_files, meta, read_ID, read_length, insert_size. 
    biqseq_norm_input_tumor = SAMTOOLS_FILTER_UNIQUE_TUMOR.out.uniq_seq.join(SAMTOOLS_STATS_INSERTSIZE_TUMOR.out.read_length_insert_size)
    // sampleID, individual_chr_seq_files, meta, read_ID, read_length, insert_size. 

    fasta_files = Channel.fromPath( file(params.ref_fa).parent + '/*_chr*' )
            .collect()
    // collect individual chr fasta files. These are located in the same directory as the main reference. 
    // if the extension of `name_chr#.fa` changes this match will break. 

    BICSEQ2_NORMALIZE_NORMAL(biqseq_norm_input_normal, fasta_files)
    BICSEQ2_NORMALIZE_TUMOR(biqseq_norm_input_tumor, fasta_files)
    // note: this can not be split by chrom, even though bicseq2 norm acts on chroms in turn, 
    // it needs all chroms to parameterize the normalization. 
    // reported error will be in these cases: "Error in bin_read: bin file is in incorrect format."

    bicseq_normal = BICSEQ2_NORMALIZE_NORMAL.out.normalized_output
        .map{it -> [it[2].patient, it[1], it[2], it[3]]}

    bicseq_tumor = BICSEQ2_NORMALIZE_TUMOR.out.normalized_output
        .map{it -> [it[2].patient, it[1], it[2], it[3]]}

    bicseq2_seg_input = bicseq_normal.cross(bicseq_tumor)
                            .map{normal, tumor -> 
                                    def meta = [:]
                                    meta.patient    = normal[2].patient
                                    meta.normal_id  = normal[2].sampleID
                                    meta.tumor_id   = tumor[2].sampleID
                                    meta.sex        = normal[2].sex
                                    meta.id         = "${tumor[2].patient}--${tumor[2].tumor_id}--${tumor[2].normal_id}".toString()

                                    ["${tumor[2].patient}--${tumor[2].tumor_id}--${tumor[2].normal_id}".toString(), normal[1], tumor[1], normal[2], normal[3], tumor[3]]}
                                    // sampleID, individual_normal_norm_bin_files, individual_tumor_norm_bin_files, metadata, norm_readID, tumor_readID. 
                                    // The metadata object here is reset following the cross. So that ID matches up again. 
                                    // It is possible that in many to many or one to many crosses, the ID field will not reflect the crossed samples. 

    BICSEQ2_SEG(bicseq2_seg_input)
    // NOTE: with insufficent coverage, the segmentation will fail because the 'lamda' factor can not be properly optimized. 

    bicseq2_tumoronly_input = BICSEQ2_NORMALIZE_TUMOR.out.normalized_output
        .filter{it[2].normal_id == 'NA12878'}

    BICSEQ2_SEG_UNPAIRED(bicseq2_tumoronly_input)
    
    bicseq2_calls = BICSEQ2_SEG_UNPAIRED.out.bicseq2_sv_calls
                    .map{it -> [it[3].id, it[1], it[2], it[3], it[4], it[5], it[6]]}
                    .mix(BICSEQ2_SEG.out.bicseq2_sv_calls)
    // remap output from unpaired bicseq2 to standard format for bicseq2 paired. And mix both channel outputs. This is passed to annotation. 

    // Step 15: MSI
    MSISENSOR2_MSI(ch_msisensor2_input)
    
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
    GRIPSS_SOMATIC_FILTER.out.gripss_filtered_bgz

    Bicseq2
    BICSEQ2_SEG.out.bicseq2_sv_calls
    */

    /*
        NOTE: 
        The call merging and annotations sections of this workflow becomes highly complex. 
        Files from each caller are passed through a set of 'merge prep' steps. 
        These steps apply various functions to manipulate the VCF header, and also calls within the VCFs. 
        Once the VCFs are prepared, a merge occurs. Following the merge, non-exonic regions are parsed out, 
        and calls in those regions are passed to Lancet for confirmation/rescue. 
        Following this, confirmed calls are used as 'support' and merged back to the full caller call set. 
        Additional manipulations are done on the VCF, and then the 'final' VCF
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
                        .groupTuple(size: 5)
                        .map{sampleID, vcf, idx, meta, normal_sample, tumor_sample, tool_list -> tuple( sampleID, vcf, idx, meta.unique()[0] )  }
                        .combine(chrom_list_noY.flatten())
    // The above collects all callers on sampleID, then maps to avoid duplication of data and to drop the tool list, which is not needed anymore. 
    // Note that this could be done using 'by' in the groupTuple statement. However, the map is still required to remove the tool list. 
    // 'size: 5' corresponds to the 5 callers used in the workflow. If additional callers are added, this must be changed. 

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
                           .combine(ch_paired_samples, by: 0)
                           .map{sampleID, bed, meta, chrom, meta2, normal_bam, normal_bai, normal_name, tumor_bam, tumor_bai, tumor_name -> tuple( sampleID, bed, meta, normal_bam, normal_bai, normal_name, tumor_bam, tumor_bai, tumor_name, chrom )  }
    // The above combines output by sampleID with BAM files. Then maps to avoid duplication of data, and set input tuples for the steps that follow.  
    // Note that "combine" here, combines each output stream from VCF_TO_BED with ch_paired_samples, keeping the scattered chrom seperate. 

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

    // 2. Add Allele Count to VCF.
    //    "Runs pileup on tumor and normal bam files to compute allele counts for bi-allelic SNV and Indel variants in VCF file and adds pileup format columns to the VCF file.""
    addAlleleCounts_confirm_input = MERGE_COLUMNS.out.mergeColumn_vcf                         
                           .combine(ch_paired_samples, by: 0)
                           .map{sampleID, vcf, meta, chrom, meta2, normal_bam, normal_bai, normal_name, tumor_bam, tumor_bai, tumor_name -> tuple( sampleID, vcf, meta, normal_bam, normal_bai, tumor_bam, tumor_bai, chrom )  }
    ADD_NYGC_ALLELE_COUNTS(addAlleleCounts_confirm_input)

    // 3. Add Final Allele Counts to VCF
    ADD_FINAL_ALLELE_COUNTS(ADD_NYGC_ALLELE_COUNTS.out.vcf)

    // 4. Filter VCF based on PON
    FILTER_PON(ADD_FINAL_ALLELE_COUNTS.out.vcf)

    // 5. Filter VCF based on gnomad and "ALL_GRCh38_sites"
    FILTER_VCF(FILTER_PON.out.vcf)

    // 6. "SnvstomnvsCountsbasedfilterAnnotatehighconf" 
    //    Parses file and converts adjacent SNVs to MNVs if they have they match the MNV_ID and called_by fields.
    SNV_TO_MNV_FINAL_FILTER(FILTER_VCF.out.vcf)

    // ** Collect and Merge Chroms. 
    num_intervals = file(params.chrom_contigs).countLines().toInteger() - 2
    // number of chrom intervals split on during the above steps. A 'value' variable used in groupTuple size statement. MT and Y are removed, hence '- 2'
    chrom_merge_input = SNV_TO_MNV_FINAL_FILTER.out.vcf
                        .groupTuple(size: num_intervals)
                        .map{sampleID, vcf, meta, chrom -> tuple( sampleID, vcf, meta.unique()[0] )  }
                        // Collect scattered chroms, remap to tuple without chrom names. 

    GATK_SORTVCF_SOMATIC(chrom_merge_input)
    REORDER_VCF_COLUMNS(GATK_SORTVCF_SOMATIC.out.vcf_idx)
    // output tuple = val(sampleID), path("*_mergePrep.vcf"), val(meta). 
    // meta = [patient:test, normal_id:test, tumor_id:test2, sex:XX, id:test2_vs_test] 
    //         This named list can be accessed in the script section prior to """ via calls like: meta.patient

    // Compress and index the merged vcf
    COMPRESS_INDEX_MERGED_VCF(REORDER_VCF_COLUMNS.out.vcf)

    // ** Annotation of somatic indels and snps

    VEP_SOMATIC(COMPRESS_INDEX_MERGED_VCF.out.compressed_vcf_tbi)
    COSMIC_ANNOTATION_SOMATIC(VEP_SOMATIC.out.vcf)
    COSMIC_CANCER_RESISTANCE_MUTATION_SOMATIC(COSMIC_ANNOTATION_SOMATIC.out.vcf)

    SNPSIFT_ANNOTATE_DBSNP_SOMATIC(COSMIC_CANCER_RESISTANCE_MUTATION_SOMATIC.out.vcf.map{it -> [it[0], it[1]]}, params.dbSNP, params.dbSNP_index, 'intermediate')
    // note: existing module requires only sampleID and VCF. input remapped to required tuple.

    somatic_finalization_input = SNPSIFT_ANNOTATE_DBSNP_SOMATIC.out.vcf.join(COSMIC_CANCER_RESISTANCE_MUTATION_SOMATIC.out.vcf).map{it -> [it[0], it[1], it[3], it[4], it[5]]}
    // re-join dbSNP ID annotated VCF output with [meta], normalID, tumorID. 

    SOMATIC_VCF_FINALIZATION(somatic_finalization_input, 'filtered')

    // ** Annotation of somatic CNV and SV

    ANNOTATE_BICSEQ2_CNV(bicseq2_calls, chrom_list_noY)
    
    // note: joining on the sampleID, metadata, tumor_name, and normal_name for
    // safety. This re-arranges the values in the channel to:
    // tuple val(sampleID), val(normal_name), val(tumor_name), file(manta_vcf), file(manta_vcf_tbi), val(meta_manta), val(manta), file(gridss_bgz), val(no_idx), val(meta_gripss), val(gridss)
    // Downstream, just including sampleID, normal_name, and tumor_name to simplify a similar join that is necessary

    merge_sv_input = MANTA.out.manta_somaticsv_tbi.join(GRIPSS_SOMATIC_FILTER.out.gripss_filtered_bgz, by : [0,4,5])
    MERGE_SV(merge_sv_input, chrom_list)
    
    ANNOTATE_SV(MERGE_SV.out.merged, "main")
    ANNOTATE_SV_SUPPLEMENTAL(MERGE_SV.out.merged_suppl, "supplemental")
    ANNOTATE_GENES_SV(ANNOTATE_SV.out.annot_sv_bedpe, "main")
    ANNOTATE_GENES_SV_SUPPLEMENTAL(ANNOTATE_SV_SUPPLEMENTAL.out.annot_sv_bedpe, "supplemental")
    
    // note: joining on the sampleID, normal_name, and tumor_name for
    // safety. This re-arranges the values in the channel to:
    // tuple val(sampleID), val(normal_name), val(tumor_name), file(bicseq_annot), file(annot_sv_genes_bedpe)

    annot_sv_cnv_input = ANNOTATE_BICSEQ2_CNV.out.bicseq_annot.join(ANNOTATE_GENES_SV.out.annot_sv_genes_bedpe, by: [0,2,3])
    ANNOTATE_SV_WITH_CNV(annot_sv_cnv_input, "main")
    
    // See notes on previous step
    annot_sv_cnv_suppl_input = ANNOTATE_BICSEQ2_CNV.out.bicseq_annot.join(ANNOTATE_GENES_SV_SUPPLEMENTAL.out.annot_sv_genes_bedpe, by: [0,2,3])
    ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL(annot_sv_cnv_suppl_input, "supplemental")
    
    FILTER_BEDPE(ANNOTATE_SV_WITH_CNV.out.sv_genes_cnv_bedpe, "main")
    FILTER_BEDPE_SUPPLEMENTAL(ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL.out.sv_genes_cnv_bedpe, "supplemental")

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(JAX_TRIMMER.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_XENOME_CLASSIFY_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK_BASERECALIBRATOR.out.table.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTWGSMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CONPAIR.out.concordance.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CONPAIR.out.contamination.collect{it[1]}.ifEmpty([]))
  
    MULTIQC (
        ch_multiqc_files.collect()
    )

}

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
            if (!sample2patient.containsKey(row.sampleID.toString())) {
                sample2patient[row.sampleID.toString()] = row.patient.toString()
            } else if (sample2patient[row.sampleID.toString()] != row.patient.toString()) {
                log.error('The sample "' + row.sampleID.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sampleID.toString()] + '" in the sample sheet.')
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
            if (!(row.patient && row.sampleID)){
                log.error "Missing field in csv file header. The csv file must have fields named 'patient' and 'sampleID'."
                System.exit(1)
            }
            [[row.patient.toString(), row.sampleID.toString()], row]
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
        if (row.sampleID)  meta.sampleID  = row.sampleID.toString()

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
            meta.id         = "${row.patient}--${row.sampleID}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)

            meta.size       = 1 // default number of splitted fastq

            return [meta.id, meta, [fastq_1, fastq_2]]

        } else {
            log.error "Missing or unknown field in csv file header. Please check your samplesheet"
            System.exit(1)
        }
    }
}