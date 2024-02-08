#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CLUMPIFY} from "${projectDir}/modules/bbmap/bbmap_clumpify"
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES}	from "${projectDir}/modules/picard/picard_markduplicates"

include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"

include {GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_NORMAL;
         GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_TUMOR} from "${projectDir}/modules/gatk/gatk_getsamplename"

include {GATK_HAPLOTYPECALLER_SV_GERMLINE} from "${projectDir}/modules/gatk/gatk_haplotypecaller_sv_germline"

include {GATK_SORTVCF_GERMLINE as GATK_SORTVCF_GERMLINE;
         GATK_SORTVCF_GERMLINE as GATK_SORTVCF_GENOTYPE} from "${projectDir}/modules/gatk/gatk_sortvcf_germline"
include {GATK_GENOTYPE_GVCF} from "${projectDir}/modules/gatk/gatk_genotype_gvcf"
include {GATK_CNNSCORE_VARIANTS} from "${projectDir}/modules/gatk/gatk_cnnscorevariants"
include {GATK_VARIANTFILTRATION_AF} from "${projectDir}/modules/gatk/gatk_variantfiltration_af"
include {BCFTOOLS_COMPRESS_INDEX} from "${projectDir}/modules/bcftools/bcftools_compress_index"
include {BCFTOOLS_SPLITMULTIALLELIC_REGIONS} from "${projectDir}/modules/bcftools/bcftools_split_multiallelic_regions"
include {VEP_GERMLINE} from "${projectDir}/modules/ensembl/varianteffectpredictor_germline_mouse"
include {BCFTOOLS_REMOVESPANNING} from "${projectDir}/modules/bcftools/bcftools_remove_spanning"

include {GATK_MUTECT2} from "${projectDir}/modules/gatk/gatk_mutect2"
include {GATK_MERGEMUTECTSTATS} from "${projectDir}/modules/gatk/gatk_mergemutectstats"
include {GATK_FILTERMUECTCALLS} from "${projectDir}/modules/gatk/gatk_filtermutectcalls"
include {LANCET} from "${projectDir}/modules/nygenome/lancet"
include {MANTA} from "${projectDir}/modules/illumina/manta"
include {STRELKA2} from "${projectDir}/modules/illumina/strelka2"
include {DELLY_CALL_SOMATIC} from "${projectDir}/modules/delly/delly_call_somatic"
include {DELLY_FILTER_SOMATIC} from "${projectDir}/modules/delly/delly_filter_somatic"
include {BCFTOOLS_BCF_TO_VCF} from "${projectDir}/modules/bcftools/bcftools_bcf_to_vcf"
include {SMOOVE_CALL} from "${projectDir}/modules/smoove/smoove_call"
include {SVABA} from "${projectDir}/modules/svaba/svaba"
include {GATK_UPDATEVCFSEQUENCEDICTIONARY as SVABA_SV_UPDATE_DICTIONARY;
         GATK_UPDATEVCFSEQUENCEDICTIONARY as SVABA_INDEL_UPDATE_DICTIONARY} from "${projectDir}/modules/gatk/gatk_updatevcfsequencedictionary"

include {DELLY_CNV_SOMATIC} from "${projectDir}/modules/delly/delly_cnv_somatic"
include {BCFTOOLS_MERGE_DELLY_CNV} from "${projectDir}/modules/bcftools/bcftools_merge_delly_cnv"
include {DELLY_CLASSIFY} from "${projectDir}/modules/delly/delly_classify"
include {BCFTOOLS_QUERY_DELLY_CNV} from "${projectDir}/modules/bcftools/bcftools_query_delly_cnv"
include {PLOT_DELLY_CNV} from "${projectDir}/modules/r/plot_delly_cnv"

include {GATK_SORTVCF as GATK_SORTVCF_MUTECT;
         GATK_SORTVCF as GATK_SORTVCF_LANCET;
         GATK_SORTVCF as GATK_SORTVCF_TOOLS;
         GATK_SORTVCF as GATK_SORTVCF_TOOLS_LANCET} from "${projectDir}/modules/gatk/gatk_sortvcf_somatic_tools"

include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP_GERMLINE;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP_SOMATIC} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"

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
include {SNV_TO_MNV_FINAL_FILTER} from "${projectDir}/modules/python/python_snv_to_mnv_final_filter"

include {GATK_SORTVCF_SOMATIC} from "${projectDir}/modules/gatk/gatk_sortvcf_somatic_merge"
include {REORDER_VCF_COLUMNS} from "${projectDir}/modules/python/python_reorder_vcf_columns"
include {COMPRESS_INDEX_MERGED_VCF} from "${projectDir}/modules/tabix/compress_merged_vcf"
include {VEP_SOMATIC} from "${projectDir}/modules/ensembl/varianteffectpredictor_somatic_mouse"
include {SOMATIC_VCF_FINALIZATION} from "${projectDir}/modules/python/python_somatic_vcf_finalization_mouse"

include {ANNOTATE_DELLY_CNV} from "${projectDir}/modules/r/annotate_delly_cnv"

include {MERGE_SV} from "${projectDir}/modules/r/merge_sv_mouse"
include {ANNOTATE_SV;
         ANNOTATE_SV as ANNOTATE_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/annotate_sv_mouse"
include {ANNOTATE_GENES_SV;
         ANNOTATE_GENES_SV as ANNOTATE_GENES_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/annotate_genes_sv_mouse"
include {ANNOTATE_SV_WITH_CNV;
         ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL} from "${projectDir}/modules/r/annotate_sv_with_cnv_mouse"
include {FILTER_BEDPE;
         FILTER_BEDPE as FILTER_BEDPE_SUPPLEMENTAL} from "${projectDir}/modules/r/filter_bedpe_mouse"

include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

workflow MM_PTA {
    take:
        concat_ch
    main:
        concat_ch.map{it -> [it[0], it[2]]}.set{read_ch}
        concat_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
        
        // Optional Step -- Clumpify
        if (params.deduplicate_reads) {
            CLUMPIFY(read_ch)
            trimmer_input = CLUMPIFY.out.clumpy_fastq
        } else {
            trimmer_input = read_ch
        }
        // ** Trimmer
        FASTP(CLUMPIFY.out.clumpy_fastq)
        
        FASTQC(FASTP.out.trimmed_fastq)
        
        // ** Get Read Group Information
        READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")

        bwa_mem_mapping = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
        
        // ** BWA-MEM Alignment
        BWA_MEM(bwa_mem_mapping)

        // ** Sort mapped reads
        PICARD_SORTSAM(BWA_MEM.out.sam)

        // ** Markduplicates
        PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

        PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai).join(meta_ch).branch{
            normal: it[3].status == 0
            tumor:  it[3].status == 1
        }.set{ch_final_bam}
        // re-join the sampleID to metadata information. Split normal and tumor samples into 2 different paths. 
        // Process tumor and normal BAMs seperately for conpair. For calling, use mapped and crossed data. 

        // ** Get alignment and WGS metrics
        PICARD_COLLECTALIGNMENTSUMMARYMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)
        PICARD_COLLECTWGSMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)


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

        // Restore un-paired tumor samples, and add the proxy normal sample as pairing in those cases
        ch_paired_samples = ch_tumor_to_cross
            .mix(ch_paired_samples)
            .map{it -> [it[1].patient, it[1], it[2], it[3], it[4]]}.groupTuple().filter{it[2].size() == 1} 
                        // it[0] = sampleID, it[1] = meta, it[2] = bam, it[3] = bai, it[4] = sampleReadID. 
                        // unknown group size, no 'size' statement can be used in groupTuple
            .map{tumor -> 
            def meta = [:]
                meta.patient    = tumor[1][0].patient
                meta.normal_id  = params.proxy_normal_sampleName
                meta.tumor_id   = tumor[1][0].sampleID
                meta.sex        = tumor[1][0].sex
                meta.id         = "${meta.patient}--${meta.tumor_id}--${meta.normal_id}".toString()
            
                [meta.id, meta, params.proxy_normal_bam, params.proxy_normal_bai, params.proxy_normal_sampleName, tumor[2][0], tumor[3][0], tumor[4][0]]
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
            To recover un-paired tumors, pair them with the proxy normal sample and pass them with paired samples downstream,
            the group of all tumor samples: 'ch_tumor_to_cross' is mixed with the paired sample: 'ch_paired_samples'.
            Cases where tumor samples were paired are then filtered. Tumors that were paired in the cross will appear > 2 times in the mix results, and are removed via the 'filter it[2].size()==1' statement. 
            The resulting tumor-only samples are mapped into the format seen in the cross statement, with the proxy normal sample being added via parameters as the 'normal' sample. 
        */


        ch_ind_samples = ch_paired_samples
            .filter{it[4] != params.proxy_normal_sampleName}
            .multiMap{it -> 
                    normal: ["${it[1].patient}--${it[1].normal_id}".toString(), it[1], it[2], it[3], it[4]]
                    tumor:  ["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]
                    }
            ch_normal_samples = ch_ind_samples.normal.unique{it[0]}
            ch_tumor_samples  = ch_ind_samples.tumor.unique{it[0]}

        ch_tumor_only = ch_paired_samples
            .filter{it[4] == params.proxy_normal_sampleName}
            .map{it -> ["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]}
            .unique{it[0]}


        /*
            The above establishes channels needed for germline calling, and delly. 
            Those steps require BAM, index and readID. 
            Here sampleID is reset to the original ID from the CSV parser which is: 'patient--sample'
            Note that the proxy normal sample is filtered for germline and can be filtered for delly.
            Note: MSIsensor2 is unsupported for mouse, as it is lacking model files.  
        */

        // ** NEXTFLOW OPERATORS SETUP END


        // ** Germline Calling and Annotation
        
        // Find the paths of all `scattered.interval_list` files, and make tuples with an index value. 
        // This is used for for HaplotypeCaller variant regions
        
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
        chrom_channel = ch_normal_samples.combine(intervals).filter{it[4] != params.proxy_normal_sampleName}

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

        // Allele frequency and call refinement filtering. 
        GATK_VARIANTFILTRATION_AF(GATK_SORTVCF_GENOTYPE.out.vcf_idx)
        BCFTOOLS_COMPRESS_INDEX(GATK_VARIANTFILTRATION_AF.out.vcf)

        // Germline annotation - Filtered
        // 1. SplitMultiAllelicRegions & compress & index
        BCFTOOLS_SPLITMULTIALLELIC_REGIONS(BCFTOOLS_COMPRESS_INDEX.out.vcf_idx, chrom_list)
        // 2. vepPublicSvnIndel
        VEP_GERMLINE(BCFTOOLS_SPLITMULTIALLELIC_REGIONS.out.vcf_idx)
        // 3. RemoveSpanning
        BCFTOOLS_REMOVESPANNING(VEP_GERMLINE.out.vcf)
        // 5. Add dbsnpIDs
        SNPSIFT_ANNOTATE_DBSNP_GERMLINE(BCFTOOLS_REMOVESPANNING.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
        
        // ** Somatic Calling

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

        // ** Lumpy via Smoove - SV Calling
        SMOOVE_CALL(ch_paired_samples)

        // ** Delly - SV Calling
        DELLY_CALL_SOMATIC(ch_paired_samples)
        DELLY_FILTER_SOMATIC(DELLY_CALL_SOMATIC.out.bcf_csi)
        BCFTOOLS_BCF_TO_VCF(DELLY_FILTER_SOMATIC.out.bcf_csi)

        // ** Delly - CNV Calling
        DELLY_CNV_SOMATIC(ch_paired_samples)
        BCFTOOLS_MERGE_DELLY_CNV(DELLY_CNV_SOMATIC.out.bcfs)
        DELLY_CLASSIFY(BCFTOOLS_MERGE_DELLY_CNV.out.merged_bcf)
        BCFTOOLS_QUERY_DELLY_CNV(DELLY_CLASSIFY.out.bcf_csi)

        r_plot_input = DELLY_CNV_SOMATIC.out.tumor_cov.join(BCFTOOLS_QUERY_DELLY_CNV.out.segmentation_file)
                       .map{it -> tuple(it[0], it[1], it[2])}

        PLOT_DELLY_CNV(r_plot_input)

        SVABA(ch_paired_samples)
        SVABA_SV_UPDATE_DICTIONARY(SVABA.out.svaba_somatic_sv_vcf_tbi, 'svaba') 
        SVABA_INDEL_UPDATE_DICTIONARY(SVABA.out.svaba_somatic_indel_vcf_tbi, 'svaba')
        // Note: SVABA jumbles the dict header in the VCF, and must be adjusted prior to use in GATK tools.

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

        Lumpy via SMOOVE_CALL
        SMOOVE_CALL.out.vcf_tbi

        Delly SV
        BCFTOOLS_BCF_TO_VCF.out.vcf_tbi

        Delly CNV
        BCFTOOLS_QUERY_DELLY_CNV.out.segmentation_file

        SVABA_SV
        SVABA_SV_UPDATE_DICTIONARY.out.vcf_tbi

        SVABA_INDEL
        SVABA_INDEL_UPDATE_DICTIONARY.out.vcf_tbi

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

        somatic_caller_concat = MANTA.out.manta_somaticsv_tbi.concat(STRELKA2.out.strelka_snv_vcf_tbi, 
                                                                    STRELKA2.out.strelka_indel_vcf_tbi, 
                                                                    GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 
                                                                    SVABA_INDEL_UPDATE_DICTIONARY.out.vcf_tbi,
                                                                    GATK_SORTVCF_LANCET.out.vcf_tbi)

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
                            .groupTuple(size: 6)
                            .map{sampleID, vcf, idx, meta, normal_sample, tumor_sample, tool_list -> tuple( sampleID, vcf, idx, meta.unique()[0] )  }
                            .combine(chrom_list.flatten())

        // The above collects all callers on sampleID, then maps to avoid duplication of data and to drop the tool list, which is not needed anymore. 
        // Note that this could be done using 'by' in the groupTuple statement. However, the map is still required to remove the tool list. 
        // 'size: 6' corresponds to the 6 callers used in the workflow. If additional callers are added, this must be changed. 

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

        // 4. Filter VCF based on PON is not possible in mouse. We lack a PON

        // 5. Filter VCF based on gnomad and "ALL_GRCh38_sites" is not possible in mouse, as these are human only reference files. 

        // 6. "SnvstomnvsCountsbasedfilterAnnotatehighconf" 
        //    Parses file and converts adjacent SNVs to MNVs if they have they match the MNV_ID and called_by fields.
        SNV_TO_MNV_FINAL_FILTER(ADD_FINAL_ALLELE_COUNTS.out.vcf)

        // ** Collect and Merge Chroms. 
        num_intervals = file(params.chrom_contigs).countLines().toInteger() - 1
        // number of chrom intervals split on during the above steps. A 'value' variable used in groupTuple size statement. MT is removed, hence '- 1'. If dynamic 'Y' is ever implimented. This needs adjustment. 
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

        SNPSIFT_ANNOTATE_DBSNP_SOMATIC(VEP_SOMATIC.out.vcf.map{it -> [it[0], it[1]]}, params.dbSNP, params.dbSNP_index, 'intermediate')
        // note: existing module requires only sampleID and VCF. input remapped to required tuple.

        somatic_finalization_input = SNPSIFT_ANNOTATE_DBSNP_SOMATIC.out.vcf.join(VEP_SOMATIC.out.vcf).map{it -> [it[0], it[1], it[3], it[4], it[5]]}
        // re-join dbSNP ID annotated VCF output with [meta], normalID, tumorID. 

        SOMATIC_VCF_FINALIZATION(somatic_finalization_input, 'filtered')


        // ** Annotation of somatic CNV and SV

        ANNOTATE_DELLY_CNV(BCFTOOLS_QUERY_DELLY_CNV.out.segmentation_file, chrom_list)


        // note: joining on the sampleID, metadata, tumor_name, and normal_name for
        // safety. This re-arranges the values in the channel to:
        // tuple val(sampleID), val(normal_name), val(tumor_name), file(manta_vcf), file(manta_vcf_tbi), val(meta_manta), val(manta), file(smoove_vcf), val(smoove_tbi), val(meta_smoove), val(lumpy) # NOTE: SMOOVE runs lumpy, svtyper and duphold. 
        // to then join delly, the values are re-mapped to place tumor name and normal name back into the 4 and 5 index positions. 
        // caller names and duplicated repeat metadata tuples are removed. 
        // Downstream, just including sampleID, normal_name, and tumor_name to simplify a similar join that is necessary

        merge_sv_input = MANTA.out.manta_somaticsv_tbi.join(SMOOVE_CALL.out.vcf_tbi, by : [0,4,5])
                         .map{it -> tuple(it[0], it[3], it[4], it[5], it[1], it[2], it[7], it[8])}

                         .join(BCFTOOLS_BCF_TO_VCF.out.vcf_tbi, by : [0,4,5])
                         .map{it -> tuple(it[0], it[3], it[4], it[6], it[1], it[2], it[7], it[8], it[9])}
                         .join(SVABA_SV_UPDATE_DICTIONARY.out.vcf_tbi, by : [0,4,5])
                         .map{it -> tuple(it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8], it[9], it[10])}
        // tuple val(sampleID), val(normal_name), val(tumor_name), path(manta_vcf), path(manta_vcf_tbi), path(smoove_vcf), path(smoove_tbi), path(delly_vcf), path(delly_tbi), path(svaba_vcf), path(svaba_tbi)

        MERGE_SV(merge_sv_input, chrom_list) //slop=1000 and min_sv_length=200
        
        ANNOTATE_SV(MERGE_SV.out.merged, "main")
        ANNOTATE_SV_SUPPLEMENTAL(MERGE_SV.out.merged_suppl, "supplemental")
        ANNOTATE_GENES_SV(ANNOTATE_SV.out.annot_sv_bedpe, "main")
        ANNOTATE_GENES_SV_SUPPLEMENTAL(ANNOTATE_SV_SUPPLEMENTAL.out.annot_sv_bedpe, "supplemental")
        
        // note: joining on the sampleID, normal_name, and tumor_name for
        // safety. This re-arranges the values in the channel to:
        // tuple val(sampleID), val(normal_name), val(tumor_name), file(delly_annot), file(annot_sv_genes_bedpe)

        annot_sv_cnv_input = ANNOTATE_DELLY_CNV.out.delly_annot.join(ANNOTATE_GENES_SV.out.annot_sv_genes_bedpe, by: [0,2,3])
        ANNOTATE_SV_WITH_CNV(annot_sv_cnv_input, "main")
        
        // See notes on previous step
        annot_sv_cnv_suppl_input = ANNOTATE_DELLY_CNV.out.delly_annot.join(ANNOTATE_GENES_SV_SUPPLEMENTAL.out.annot_sv_genes_bedpe, by: [0,2,3])
        ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL(annot_sv_cnv_suppl_input, "supplemental")
        
        FILTER_BEDPE(ANNOTATE_SV_WITH_CNV.out.sv_genes_cnv_bedpe, "main")
        FILTER_BEDPE_SUPPLEMENTAL(ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL.out.sv_genes_cnv_bedpe, "supplemental")
   
        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTWGSMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
    
        MULTIQC (
            ch_multiqc_files.collect()
        )

}
