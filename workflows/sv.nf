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
include {GATK_SORTVCF as GATK_SORTVCF_GERMLINE;
         GATK_SORTVCF as GATK_SORTVCF_GENOTYPE;
         GATK_SORTVCF as GATK_SORTVCF_MUTECT;
         GATK_SORTVCF as GATK_SORTVCF_LANCET} from "${projectDir}/modules/gatk/gatk_sortvcf"
include {GATK_GENOTYPE_GVCF} from "${projectDir}/modules/gatk/gatk_genotype_gvcf"
include {GATK_CNNSCORE_VARIANTS} from "${projectDir}/modules/gatk/gatk_cnnscorevariants"
include {GATK_FILTER_VARIANT_TRANCHES} from "${projectDir}/modules/gatk/gatk_filtervarianttranches"
include {GATK_VARIANTFILTRATION_AF} from "${projectDir}/modules/gatk/gatk_variantfiltration_af"
include {BCFTOOLS_GERMLINE_FILTER} from "${projectDir}/modules/bcftools/bcftools_germline_filter"
include {GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_NORMAL;
         GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_TUMOR} from "${projectDir}/modules/gatk/gatk_getsamplename"
include {GATK_MUTECT2} from "${projectDir}/modules/gatk/gatk_mutect2"
include {GATK_MERGEMUTECTSTATS} from "${projectDir}/modules/gatk/gatk_mergemutectstats"
include {GATK_FILTERMUECTCALLS} from "${projectDir}/modules/gatk/gatk_filtermutectcalls"
include {MANTA} from "${projectDir}/modules/illumina/manta"
include {STRELKA2} from "${projectDir}/modules/illumina/strelka2"
include {LANCET} from "${projectDir}/modules/nygenome/lancet"
include {GRIDSS_PREPROCESS} from "${projectDir}/modules/gridss/gridss_preprocess"
include {GRIDSS_ASSEMBLE} from "${projectDir}/modules/gridss/gridss_assemble"
include {GRIDSS_CALLING} from "${projectDir}/modules/gridss/gridss_calling"
include {GRIDSS_CHROM_FILTER} from "${projectDir}/modules/gridss/gridss_chrom_filter"
include {GRIDSS_SOMATIC_FILTER} from "${projectDir}/modules/gridss/gridss_somatic_filter"
include {SAMTOOLS_STATS_INSERTSIZE as SAMTOOLS_STATS_INSERTSIZE_NORMAL;
         SAMTOOLS_STATS_INSERTSIZE as SAMTOOLS_STATS_INSERTSIZE_TUMOR} from "${projectDir}/modules/samtools/samtools_stats_insertsize"
include {SAMTOOLS_FILTER_UNIQUE as SAMTOOLS_FILTER_UNQIUE_NORMAL;
         SAMTOOLS_FILTER_UNIQUE as SAMTOOLS_FILTER_UNQIUE_TUMOR} from "${projectDir}/modules/samtools/samtools_filter_unique_reads"
include {BICSEQ2_NORMALIZE as BICSEQ2_NORMALIZE_NORMAL;
         BICSEQ2_NORMALIZE as BICSEQ2_NORMALIZE_TUMOR} from "${projectDir}/modules/biqseq2/bicseq2_normalize"
include {BICSEQ2_SEG} from "${projectDir}/modules/biqseq2/bicseq2_seg"
include {SVABA} from "${projectDir}/modules/svaba/svaba"
include {LUMPY_SV} from "${projectDir}/modules/lumpy_sv/lumpy_sv"

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

    // Step 13: Germline Calling
    
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


    // NEED TO ADD ANNOTATION OF GERMLINE. 


    // Step 14: Somatic Calling

    // Read a list of contigs from parameters to provide to GATK as intervals
    chroms = Channel
        .fromPath("${params.chrom_contigs}")
        .splitText()
        .map{it -> it.trim()}

    // Get a list of primary chromosomes and exclude chrM (dropRight(1))
    chrom_list = chroms.collect().dropRight(1)

    // Applies scatter intervals from above to the BQSR bam file
    somatic_calling_channel = ch_cram_variant_calling_pair.combine(intervals)

    // // Applies scatter intervals from above to the BQSR bam file
    // somatic_calling_channel = ch_cram_variant_calling_pair.combine(chroms)
    // // NOTE: The above code line will split by Mutect2 calling by indivdiaul chromosomes 'chr'. 
    // //       Entire chromosomes are scattered. For WGS, this is computationally intensive. 
    // //       We changed to calling to be done based on the same intervals passed to the germline caller. 
    // //       If complete chromosomes are requried, the above line of code can be uncommented. 

    // Mutect2
    // STEPS: Call on each chromosome / interval. 
    //        Prior to 'filtermutectcalls' vcfs must be merged. (GATK best practice) 
    //        Prior to 'filtermutectcalls' "stats" files from mutect2 must be merged. (GATK best practice) 
    //        Merge vcfs and stats must be Nextflow joined prior to 'filtermutectcalls' to avoid samples being confounded. 

    GATK_MUTECT2(somatic_calling_channel)
    GATK_SORTVCF_MUTECT(GATK_MUTECT2.out.vcf.groupTuple(), 'vcf')
    GATK_MERGEMUTECTSTATS(GATK_MUTECT2.out.stats.groupTuple())
    filter_mutect_input = GATK_SORTVCF_MUTECT.out.vcf_idx.join(GATK_MERGEMUTECTSTATS.out.stats)
    GATK_FILTERMUECTCALLS(filter_mutect_input)
    // additional NYGC steps not used: add commands to VCF, and reorder VCF columns. 

    // Manta
    MANTA(ch_cram_variant_calling_pair)
    // additional NYGC steps not used: add commands to VCF, and reorder VCF columns. 
    // FilterNonpass is used in NYGC with `SelectVariants` and `--exclude-filtered`. Do we want hard filtering? 

    // Strelka2
    strekla2_input = ch_cram_variant_calling_pair.join(MANTA.out.manta_smallindel_vcf_tbi)
    STRELKA2(strekla2_input)
    // additional NYGC steps not used: add commands to VCF, and reorder VCF columns. 

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
    // additional NYGC steps not used: add commands to VCF, and reorder VCF columns. 
    GATK_SORTVCF_LANCET(LANCET.out.lancet_vcf.groupTuple(), 'vcf')
    //NOTE: :: ::: ::: :: We need to be able to capture this output. 

    // Gridss
    GRIDSS_PREPROCESS(ch_cram_variant_calling_pair)
    gridss_assemble_input = ch_cram_variant_calling_pair.join(GRIDSS_PREPROCESS.out.gridss_preproc)
    GRIDSS_ASSEMBLE(gridss_assemble_input)
    gridss_call_input = ch_cram_variant_calling_pair.join(GRIDSS_ASSEMBLE.out.gridss_assembly)
    GRIDSS_CALLING(gridss_call_input)
    GRIDSS_CHROM_FILTER(GRIDSS_CALLING.out.gridss_vcf, chrom_list)
    GRIDSS_SOMATIC_FILTER(GRIDSS_CHROM_FILTER.out.gridss_chrom_vcf, params.gridss_pon)
    // gridss somatic filter will need a higher coverage dataset for testing. 
    // additional NYGC steps not used: add commands to VCF, and reorder VCF columns. 

    // BicSeq2
    SAMTOOLS_STATS_INSERTSIZE_NORMAL(ch_bam_normal_to_cross)
    SAMTOOLS_STATS_INSERTSIZE_TUMOR(ch_bam_tumor_to_cross)

    SAMTOOLS_FILTER_UNQIUE_NORMAL(ch_bam_normal_to_cross, chrom_list)
    SAMTOOLS_FILTER_UNQIUE_TUMOR(ch_bam_tumor_to_cross, chrom_list)

    biqseq_norm_input_normal = SAMTOOLS_FILTER_UNQIUE_NORMAL.out.uniq_seq.join(SAMTOOLS_STATS_INSERTSIZE_NORMAL.out.read_length_insert_size)
    // sampleID, individual_chr_seq_files, read_ID, read_length, insert_size. 
    biqseq_norm_input_tumor = SAMTOOLS_FILTER_UNQIUE_TUMOR.out.uniq_seq.join(SAMTOOLS_STATS_INSERTSIZE_TUMOR.out.read_length_insert_size)
    // sampleID, individual_chr_seq_files, read_ID, read_length, insert_size. 

    fasta_files = Channel.fromPath( file(params.ref_fa).parent + '/*_chr*' )
            .collect()
    // collect individual chr fasta files. These are located in the same directory as the main reference. 
    // if the extension of `name_chr#.fa` changes this match will break. 
    // Can this be made more flexible without requiring another input parameter?

    BICSEQ2_NORMALIZE_NORMAL(biqseq_norm_input_normal, fasta_files)
    BICSEQ2_NORMALIZE_TUMOR(biqseq_norm_input_tumor, fasta_files)
    // bicseq2 normalize will need a higher coverage dataset for testing. 

    bicseq2_seg_input = BICSEQ2_NORMALIZE_NORMAL.out.normalized_output.join(BICSEQ2_NORMALIZE_TUMOR.out.normalized_output)
    // sampleID, individual_normal_norm_bin_files, individual_tumor_norm_bin_files. 
    BICSEQ2_SEG(bicseq2_seg_input)

    // Svaba
    SVABA(ch_cram_variant_calling_pair)

    // Lumpy
    LUMPY_SV(ch_cram_variant_calling_pair)

    // Step NN: Get alignment and WGS metrics
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(GATK_APPLYBQSR.out.bam)
    PICARD_COLLECTWGSMETRICS(GATK_APPLYBQSR.out.bam)

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