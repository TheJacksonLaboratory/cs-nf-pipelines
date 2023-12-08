#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/somatic_wes.nf"
include {param_log} from "${projectDir}/bin/log/somatic_wes.nf"

include {CONCATENATE_PTA_FASTQ} from "${projectDir}/subworkflows/concatenate_pta_fastq"

include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {XENOME_CLASSIFY} from "${projectDir}/modules/xenome/xenome"
// include {GZIP} from "${projectDir}/modules/utility_modules/gzip"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"

include {ANCESTRY} from "${projectDir}/workflows/ancestry"

include {GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_NORMAL;
         GATK_GETSAMPLENAME as GATK_GETSAMPLENAME_TUMOR} from "${projectDir}/modules/gatk/gatk_getsamplename"

include {CNV} from "${projectDir}/subworkflows/somatic_wes_pta_cnv"

include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration_mutect2"
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_MUTECT2} from "${projectDir}/modules/gatk/gatk_mutect2_wes_pta"

include {GATK_GETPILEUPSUMMARIES as GATK_GETPILEUPSUMMARIES_NORMAL;
         GATK_GETPILEUPSUMMARIES as GATK_GETPILEUPSUMMARIES_TUMOR} from "${projectDir}/modules/gatk/gatk_getpileupsummaries"
include {GATK_CALCULATECONTAMINATION} from "${projectDir}/modules/gatk/gatk_calculatecontamination"

include {GATK_LEARNREADORIENTATIONMODEL} from "${projectDir}/modules/gatk/gatk_learnreadorientationmodel"

include {GATK_FILTERMUECTCALLS} from "${projectDir}/modules/gatk/gatk_filtermutectcalls_wes"

include {MSISENSOR2_MSI} from "${projectDir}/modules/msisensor2/msisensor2"

include {GATK_MERGEVCF as GATK_MERGEVCF_UNANNOTATED;
         GATK_MERGEVCF as GATK_MERGEVCF_ANNOTATED} from "${projectDir}/modules/gatk/gatk_mergevcf"

include {BEDOPS_SORT} from "${projectDir}/modules/bedops/bedops_sort"
include {BEDOPS_WINDOW} from "${projectDir}/modules/bedops/bedops_window"
include {TMB_SCORE} from "${projectDir}/modules/tumor_mutation_burden/tmb_score"

include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from "${projectDir}/modules/cosmic/cosmic_annotation"
include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_DBSNP} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {PICARD_COLLECTHSMETRICS} from "${projectDir}/modules/picard/picard_collecthsmetrics"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel


if (!params.csv_input) {
    exit 1, "CSV Input Required"
}

if (params.gen_org == 'mouse' && params.pdx) {
    exit 1, "PDX workflow was called; however, `--gen_org` was set to: ${params.gen_org}. This is an invalid parameter combination. `--gen_org` must == 'human' for PDX analysis."
}

if (params.gen_org == 'mouse') {
    exit 1, "`--gen_org` was set to: ${params.gen_org}. Somatic WES currently supports only human data. `--gen_org` must == 'human' for this analysis."
}


workflow SOMATIC_WES_PTA {

    if (params.csv_input) {
        ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
        // Concat local Fastq files from CSV input if required.
            CONCATENATE_PTA_FASTQ(ch_input_sample)
            
    }

    CONCATENATE_PTA_FASTQ.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
    CONCATENATE_PTA_FASTQ.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}

    // Step 1: Read Trim
    FASTP(read_ch)

    xenome_input = FASTP.out.trimmed_fastq
    
    FASTQC(FASTP.out.trimmed_fastq)

    // Step 3: Get Read Group Information
    READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")

    // Step 1a: Xenome if PDX data used.
    ch_XENOME_CLASSIFY_multiqc = Channel.empty() //optional log file. 
    if (params.pdx){
        FASTP.out.trimmed_fastq.join(meta_ch).branch{
            normal: it[2].status == 0
            tumor:  it[2].status == 1
        }.set{fastq_files}

        normal_fastqs = fastq_files.normal.map{it -> [it[0], it[1]] }

        // Xenome Classification
        XENOME_CLASSIFY(fastq_files.tumor.map{it -> [it[0], it[1]] })
        ch_XENOME_CLASSIFY_multiqc = XENOME_CLASSIFY.out.xenome_stats //set log file for multiqc

        // GZIP(XENOME_CLASSIFY.out.xenome_human_fastq)

        // Step 4: BWA-MEM Alignment
        bwa_mem_mapping = XENOME_CLASSIFY.out.xenome_human_fastq.mix(normal_fastqs).join(READ_GROUPS.out.read_groups)

    } else { 
        bwa_mem_mapping = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
    }

    BWA_MEM(bwa_mem_mapping)

    // Step 5: Variant Preprocessing - Part 1
    PICARD_SORTSAM(BWA_MEM.out.sam)
    PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

    // Step 6: Variant Pre-Processing - Part 2
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)

    apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)
    GATK_APPLYBQSR(apply_bqsr)

    // Step 7: QC Metrics
    collect_metrics = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
    PICARD_COLLECTHSMETRICS(collect_metrics)

    // Nextflow channel processing
    // https://github.com/nf-core/sarek/blob/master/workflows/sarek.nf#L854

    GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai).join(meta_ch).branch{
        normal: it[3].status == 0
        tumor:  it[3].status == 1
    }.set{ch_final_bam}
    // re-join the sampleID to metadata information. Split normal and tumor samples into 2 different paths. 
    // Process tumor and normal BAMs seperately as needed. For calling, use mapped and crossed data. 

    ANCESTRY(ch_final_bam.normal.map{ it -> [it[0], it[1], it[2]]})

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

    ch_ind_samples = ch_paired_samples
        .multiMap{it -> 
                normal: ["${it[1].patient}--${it[1].normal_id}".toString(), it[1], it[2], it[3], it[4]]
                tumor:  ["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]
                }
        ch_normal_samples = ch_ind_samples.normal.unique{it[0]}
        ch_tumor_samples  = ch_ind_samples.tumor.unique{it[0]}

    ch_msisensor2_input = ch_paired_samples
        .map{["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6], it[7]]}
        .unique{it[0]}

    // Step: MSI
    MSISENSOR2_MSI(ch_msisensor2_input)

    // Step: CNV and HRD
    CNV(ch_paired_samples)

    // Sample Contamination Analysis
    GATK_GETPILEUPSUMMARIES_NORMAL(ch_normal_samples)
    GATK_GETPILEUPSUMMARIES_TUMOR(ch_tumor_samples)

    contam_input = GATK_GETPILEUPSUMMARIES_NORMAL.out.pileup_summary
                    .map{it -> [it[1].patient, it[1], it[2]]}
                    .cross(
                        tumor_pileups = GATK_GETPILEUPSUMMARIES_TUMOR.out.pileup_summary
                        .map{it -> [it[1].patient, it[1], it[2]]}
                    )
                    .map{normal, tumor -> [tumor[1].id, normal[2], tumor[2]]}

    GATK_CALCULATECONTAMINATION(contam_input)

    // ** Variant Calling

    // // Step: Mutect2
    GATK_MUTECT2(ch_paired_samples)

    if (params.ffpe) {
        GATK_LEARNREADORIENTATIONMODEL(GATK_MUTECT2.out.f1r2)
        filtermutectcalls_input = GATK_MUTECT2.out.vcf_tbi_stats.join(GATK_CALCULATECONTAMINATION.out.contam_segments).join(GATK_LEARNREADORIENTATIONMODEL.out.model_file)
    } else {
        filtermutectcalls_input = GATK_MUTECT2.out.vcf_tbi_stats.join(GATK_CALCULATECONTAMINATION.out.contam_segments)
                                  .map{it -> [it[0], it[1], it[2], it[3], it[4], it[5], 'no_read_model']}
    }

    GATK_FILTERMUECTCALLS(filtermutectcalls_input)

    // Step 8: Variant Filtration
    // SNP
    GATK_SELECTVARIANTS_SNP(GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 'SNP', 'selected_SNP')

    var_filter_snp = GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.idx)
    GATK_VARIANTFILTRATION_SNP(var_filter_snp, 'SNP')

    // INDEL
    GATK_SELECTVARIANTS_INDEL(GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 'INDEL', 'selected_INDEL')

    var_filter_indel = GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.idx)
    GATK_VARIANTFILTRATION_INDEL(var_filter_indel, 'INDEL')

    // Step 9: Post Variant Calling Processing - Part 1
    // 
    SNPSIFT_ANNOTATE_SNP_DBSNP(GATK_VARIANTFILTRATION_SNP.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
    SNPSIFT_ANNOTATE_SNP_COSMIC(SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
    SNPEFF_SNP(SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf, 'SNP', 'vcf')
    SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
    SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')

    // INDEL
    SNPSIFT_ANNOTATE_INDEL_DBSNP(GATK_VARIANTFILTRATION_INDEL.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
    SNPSIFT_ANNOTATE_INDEL_COSMIC(SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
    SNPEFF_INDEL(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf, 'INDEL', 'vcf')
    SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
    SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')

    // Step 10: Post Variant Calling Processing - Part 2
    vcf_files_unannotated = SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf.join(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf)
    GATK_MERGEVCF_UNANNOTATED(vcf_files_unannotated, 'SNP_INDEL_filtered_unannotated_final')

    BEDOPS_SORT(params.target_gatk)
    BEDOPS_WINDOW(BEDOPS_SORT.out.sorted_bed, params.hg38_windows)

    tmb_input = GATK_MERGEVCF_UNANNOTATED.out.vcf.join(ch_paired_samples)
                .map{it -> [it[0], it[1], "${it[2].patient}--${it[2].tumor_id}"]}

    TMB_SCORE(tmb_input, BEDOPS_WINDOW.out.window_bed)

    vcf_files_annotated = SNPEFF_ONEPERLINE_SNP.out.vcf.join(SNPEFF_ONEPERLINE_INDEL.out.vcf)
    GATK_MERGEVCF_ANNOTATED(vcf_files_annotated, 'SNP_INDEL_filtered_annotated_final')
    
    SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF_ANNOTATED.out.vcf)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK_BASERECALIBRATOR.out.table.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTHSMETRICS.out.hsmetrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_XENOME_CLASSIFY_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK_FILTERMUECTCALLS.out.stats.collect{it[1]}.ifEmpty([]))

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
