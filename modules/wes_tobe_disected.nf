#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show()


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Configurable variable parameters specific to individual runs:
    params.fastqInputs          = null // absolute path to input fastq(s)
    params.outdir               = null // directory where output will be stored
    params.tmpdir               = null // directory where files will be stored while pipeline is running
    params.ref_fa               = null // absolute path to reference genome .fa file
    params.target_gatk          = null // absolute path to gatk target bed file
    params.target_picard        = null // absolute path to picard target bed file
    params.dbSNP                = null // absolute path to dbSNP file
    params.snpEff_config        = null // absolute path to snpEff config file
}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow -c path/to/params.config run path/to/pipeline.nf -profile slurm, singularity
            (The params.config file needs to have the following mandatory parameters
             OR they need to specified on the command line.)

    Mandatory:
        --fastqInputs             absolute path to input fastq(s).
        --outdir                  directory where output will be stored.
        --tmpdir                  directory where temporary files will be held.
        --ref_fa                  absolute path to reference genome .fa file.
        --target_gatk             absolute path to gatk target bed file.
        --target_picard           absolute path to picard target bed file.
        --dbSNP                   absolute path to dbSNP file.
        --snpEff_config           absolute path to snpEff config file.

    Optional:
        -name                   Name for the pipeline run. If not specified Nextflow will
                                automatically generate a random mnemonic.

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
////// Required parameters \\\\\\
if ( ! params.fastqInputs ) {
    exit 1, "Parameter ERROR: fastqInputs ($params.fastqInputs) must be the absolute path to fastq file(s) with _R1 and _R2 fastq filenames."
}
if ( ! params.outdir ) {
    exit 1, "Parameter ERROR: Output directory parameter must be specified."
}
if ( ! file(params.outdir).exists() ) {
    exit 1, "Parameter ERROR: Missing output directory ($params.outdir) is not found: check if path is correct."
}
if ( ! params.tmpdir ) {
    exit 1, "Parameter ERROR: Temporary directory parameter must be specified."
}
if ( ! file(params.tmpdir).exists() ) {
    exit 1, "Parameter ERROR: Missing temporary directory ($params.tmpdir) is not found: check if path is correct."
}
if ( ! params.ref_fa ) {
    exit 1, "Parameter ERROR: absolute path to reference genome .fa file must be specified."
}
if ( ! params.target_gatk ) {
    exit 1, "Parameter ERROR: absolute path to gatk target bed file must be specified."
}
if ( ! params.target_picard ) {
    exit 1, "Parameter ERROR: absolute path to picard target bed file must be specified."
}
if ( ! params.dbSNP ) {
    exit 1, "Parameter ERROR: absolute path to dbSNP file must be specified."
}
if ( ! params.snpEff_config ) {
    exit 1, "Parameter ERROR: absolute path to snpEff config file must be specified."
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
println "Pipeline         = ${workflow.manifest.name}"
println "Description      = ${workflow.manifest.description}"
if(workflow.revision) {
    println "Pipeline Release = ${workflow.revision}"
}
println "Run Name         = ${workflow.runName}"
println "User             = ${workflow.userName}"
println "Config Profile   = ${workflow.profile}"
println "Config Files     = ${workflow.configFiles}"
println "Command Line     = ${workflow.commandLine}"
println "Nextflow Info    = v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
println "Launch dir       = ${workflow.launchDir}"
println "Working dir      = ${workflow.workDir}"
println "Workflow dir     = ${workflow.projectDir}"

// Pipeline Params:
println "Parameters......"                           
println ".  Output dir                              = ${params.outdir}"
println ".  Temporary dir                           = ${params.tmpdir}"
println ".  FASTQ inputs                            = ${params.fastqInputs}"
println ".  Reference genome                        = ${params.ref_fa}"
println ".  GATK target bed file                    = ${params.target_gatk}"
println ".  Picard target bed file                  = ${params.target_picard}"
println ".  dbSNP file                              = ${params.dbSNP}"
println ".  snpEff config file                      = ${params.snpEff_config}"
println "Run Start Time                             = ${workflow.start}"


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def timestamp = new Date().format("yyMMdd-HHmmss")
def tmpdir    = params.tmpdir


def fqPairRE = ~/_R1.*\..*f[ast]*q.*$/

fqin = params.fastqInputs.tokenize(",")
fqR1 = file(fqin[0])


// Need to specify absolute path to fastq file(s) because relative paths cause symlink breakage
fqR1p = fqR1.toAbsolutePath().toString()
fqR2p = file(fqin[1]).toAbsolutePath().toString()

def sampleID = ( fqR1.name - fqPairRE )


//~~~~~~~~~~ Initial Channel of SampleID and Fastq Data ~~~~
Channel.of( sampleID, fqR1p, fqR2p )
       .toList()
       .set { sample_fastqs_ch }


//~~~~~~~~~~~~~~~~ Primary publishDir == sample_tmpdir ~~~~~
def sample_tmpdir = "${params.tmpdir}/${timestamp}_${sampleID}"

// Step 1: Qual_Stat
process qual_stat {
  tag "sampleID"
  label 'long_mem'
  label 'python'

  publishDir "${sample_tmpdir}_tmp", pattern: "*fastq.gz_stat", mode: 'copy'

  input:
  tuple sampleID, read1, read2 from sample_fastqs_ch

  output:
  tuple sampleID, file("${sampleID}_R{1,2}*filtered_trimmed") into trimmed_fastq
  file "*.fastq.gz_stat"
  tuple sampleID, file("*.fastq.gz_stat") into fq_stats, dummy_fq_stats

  script:
  log.info "-----Qual_Stat running on ${sampleID}-----"

  if (params.reads == "SE"){
    mode_HQ="-S -M"
    inputfq="${read1}"
  }
  if (params.reads == "PE"){
    mode_HQ="-M"
    inputfq="${read1} ${read2}"
  }
  """

  python ${params.filter_trim} $mode_HQ ${params.min_pct_hq_reads}  $inputfq

  """
}

// Step 2: Get Read Group Information and BWA-MEM Alignment
process read_grp_bwa_mem {
  tag "sampleID"
  label 'long_high_mem'
  label 'bwa'

  input:
  tuple sampleID, file(trimmed) from trimmed_fastq

  output:
  tuple sampleID, file("*.sam") into bwa_mem_out, dummy_bwa_mem_out

  script:
  log.info "-----bwa-mem alignment running on ${sampleID}-----"

  if (params.reads == "SE"){
    inputfq="${trimmed[0]}"
    }
  if (params.reads == "PE"){
    inputfq="${trimmed[0]} ${trimmed[1]}"
    }

  """

  python ${params.read_grp_det} ${trimmed[0]} -o ${sampleID}_read_group.txt

  rg=\$(cat ${sampleID}_read_group.txt)

  /bwa-0.7.9a/bwa mem -M ${params.mismatch_penalty} -t ${task.cpus} -R \${rg} ${params.ref_fa} $inputfq > ${sampleID}.sam

  """
}

// Step 3: Variant Preprocessing 1
process variant_preproc_1 {
  tag "sampleID"
  label 'long_high_mem'
  label 'picard_2_8_1'

  publishDir "${sample_tmpdir}_tmp", pattern: "*dup_metrics.dat", mode: 'copy'

  input:
  tuple sampleID, file(sam) from bwa_mem_out

  output:
  tuple sampleID, file("*_dedup.bam") into bam_dedup
  tuple sampleID, file("*_dedup.bai") into bai_dedup
  tuple sampleID, file("*metrics.dat") into picard_metrics, dummy_picard_metrics
  file "*dup_metrics.dat"

  script:
  log.info "-----Variant pre-processing part 1 running on ${sampleID}-----"

  """

  java -Djava.io.tmpdir=$TMPDIR -Xmx24g -jar /picard.jar SortSam \
  SO=coordinate \
  INPUT=${sam} \
  OUTPUT=${sampleID}_aln.bam  \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx24g -jar /picard.jar MarkDuplicates \
  I=${sampleID}_aln.bam \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.dat \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT

  """
}

// Step 4: Variant Pre-Processing 2
process variant_preproc_2 {
  tag "sampleID"
  label 'long_high_mem'
  label 'gatk4'


  publishDir "${sample_tmpdir}_tmp", pattern: "*realigned_BQSR*", mode: 'copy'

  input:
  tuple sampleID, file(dedup_bam) from bam_dedup
  tuple sampleID, file(dedup_bai) from bai_dedup

  output:
  tuple sampleID, file("*realigned_BQSR.bam") into realigned_bam1, realigned_bam2, realigned_bam3
  tuple sampleID, file("*realigned_BQSR*bai") into realigned_bai1, realigned_bai2, realigned_bai3

  script:
  log.info "-----Variant pre-processing part 2 running on ${sampleID}-----"

  if (params.gen_org == "human") {

    """

    gatk --java-options "-Xmx24g" BaseRecalibrator \
    -I ${dedup_bam} \
    -R ${params.ref_fa} \
    --known-sites ${params.dbSNP} \
    --known-sites ${params.gold_std_indels} \
    --known-sites ${params.phase1_1000G} \
    -O recal_data.table \


    gatk --java-options "-Xmx24g" ApplyBQSR \
     -R ${params.ref_fa} \
     -I ${dedup_bam} \
     --bqsr-recal-file recal_data.table \
     -O ${sampleID}_realigned_BQSR.bam


    samtools index ${sampleID}_realigned_BQSR.bam

    rm -rf *_dedup.bam *_dedup.bai recal_data.table 

    """
    }

  else if (params.gen_org == "mouse") {

    """

    mv "${dedup_bam}" "${sampleID}_realigned_BQSR.bam"

    samtools index ${sampleID}_realigned_BQSR.bam

    rm -rf *_dedup.bam *_dedup.bai

    """
    }
}

// Step 5: Variant Preprocessing 3
process variant_preproc_3 {
  tag "sampleID"
  label 'short_med_mem'
  label 'picard_1_95'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(bam_realigned) from realigned_bam1
  tuple sampleID, file(bai_realigned) from realigned_bai1

  output:
  tuple sampleID, file("*Metrics.txt") into covmet, dummy_covmet
  file("*CoverageMetrics*")

  script:
  log.info "-----Variant pre-processing part 3 running on ${sampleID}-----"

  """

  java -Djava.io.tmpdir=$TMPDIR -jar -Xmx4g /picard-tools-1.95/CalculateHsMetrics.jar \
  TARGET_INTERVALS=${params.target_picard} \
  BAIT_INTERVALS=${params.bait_picard} \
  REFERENCE_SEQUENCE=${params.ref_fa} \
  INPUT=${bam_realigned} \
  OUTPUT=${sampleID}_CoverageMetrics.txt \
  VALIDATION_STRINGENCY=SILENT

  """
}

// Step 6a: Variant Calling
process variant_calling {
  tag "sampleID"
  label 'variant_calling_mem'
  label 'gatk4'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(bam_realigned) from realigned_bam2
  tuple sampleID, file(bai_realigned) from realigned_bai2

  output:
  tuple sampleID, file("*variants_raw.vcf") into raw_variants, dummy_raw_variants

  script:
  log.info "-----Variant calling running on ${sampleID}-----"

  """

  gatk --java-options "-Xmx12g" HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam_realigned} \
  --dbsnp ${params.dbSNP} \
  -O ${sampleID}_variants_raw.vcf \
  -L ${params.target_gatk} \
  -stand-call-conf ${params.call_val} \
  ${params.ploidy_val}

  """
}



// Step 6b: Variant Calling: Gvcf
process variant_calling_gvcf {
  tag "sampleID"
  label 'variant_calling_mem'
  label 'gatk4'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(bam_realigned) from realigned_bam3
  tuple sampleID, file(bai_realigned) from realigned_bai3

  output:
  file("*variants_raw.g.vcf.gz*")
  tuple sampleID, file("*raw.g.vcf.gz") into dummy_gvcf

  script:
  log.info "-----Variant calling running on ${sampleID}-----"

  """

  gatk --java-options "-Xmx12g" HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam_realigned} \
  -O ${sampleID}_variants_raw.g.vcf.gz \
  -ERC GVCF \
  -L ${params.target_gatk} \
  -stand-call-conf ${params.call_val} \
  ${params.ploidy_val}

  """
}

// Step 7: Variant Filtration
process variant_filtration {
  tag "sampleID"
  label 'short_med_mem'
  label 'gatk4'


  input:
  tuple sampleID, file(merged_raw_var) from raw_variants

  output:
  tuple sampleID, file("*only_snps_filtered.vcf"), file("*only_indels_filtered.vcf") into filt_var, dummy_filt_var

  script:
  log.info "-----Post variant calling processing, part 1 running on ${sampleID}-----"

  if (params.gen_org == "human") {

    """

    # Selecting SNPs from raw variants vcf
    gatk --java-options "-Xmx8g" SelectVariants \
    -R ${params.ref_fa} \
    -V ${merged_raw_var} \
    -select-type SNP \
    -O ${sampleID}_only_snps.vcf

    # Indexing SNP vcf
    gatk IndexFeatureFile \
    -I ${sampleID}_only_snps.vcf

    # Filtering SNPs
    gatk --java-options "-Xmx8g" VariantFiltration \
    -R ${params.ref_fa} \
    -V ${sampleID}_only_snps.vcf \
    -O ${sampleID}_only_snps_filtered.vcf \
    --cluster-window-size 10 \
    --filter-expression "DP < 25" --filter-name "LowCoverage" \
    --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
    --filter-expression "QUAL > 30.0 && QUAL < 50.0" --filter-name "LowQual" \
    --filter-expression "QD < 1.5" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "StrandBias"


    # Selecting Indels from raw variants vcf
    gatk --java-options "-Xmx8g" SelectVariants \
    -R ${params.ref_fa} \
    -V ${merged_raw_var} \
    -select-type INDEL \
    -O ${sampleID}_only_indels.vcf

    # Indexing indel vcf
    gatk IndexFeatureFile \
    -I ${sampleID}_only_indels.vcf


    # Filtering Indels
    gatk --java-options "-Xmx8g" VariantFiltration \
    -R ${params.ref_fa} \
    -V ${sampleID}_only_indels.vcf \
    -O ${sampleID}_only_indels_filtered.vcf \
    --cluster-window-size 10 \
    --filter-expression "DP < 25" --filter-name "LowCoverage" \
    --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
    --filter-expression "QUAL > 30.0 && QUAL < 50.0" --filter-name "LowQual" \
    --filter-expression "QD < 1.5" --filter-name "LowQD" \
    --filter-expression "FS > 200.0" --filter-name "StrandBias"

    """
    }

  else if (params.gen_org == "mouse") {

    """

    # Indexing vcf
    gatk IndexFeatureFile \
    -I ${merged_raw_var}


    gatk --java-options "-Xmx8g" VariantFiltration \
    -R ${params.ref_fa} \
    -V ${merged_raw_var} \
    -O ${sampleID}_only_snps_filtered.vcf \
    --cluster-window-size 10 \
    --filter-expression "DP < 25" --filter-name "LowCoverage" \
    --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
    --filter-expression "QUAL > 30.0 && QUAL < 50.0" --filter-name "LowQual" \
    --filter-expression "QD < 1.5" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "StrandBias"

    cat ${sampleID}_only_snps_filtered.vcf > ${sampleID}_only_indels_filtered.vcf

    """
    }

}



// Step 8: Post Variant Calling Processing, part 1
process post_var_proc1 {
  tag "sampleID"
  label 'med_high_mem'
  label 'gatk4'

  publishDir "${sample_tmpdir}_tmp", pattern: "*full*anno*", mode: 'copy'

  input:
  tuple sampleID, file(snps_filt), file(indels_filt) from filt_var

  output:
  tuple sampleID, file("*full_anno_snp*vcf") into annot_snp, dummy_annot_snp
  tuple sampleID, file("*_full_anno_indel*vcf") into annot_indel, dummy_annot_indel

  shell:
  log.info "-----Post variant calling processing, part 1 running on ${sampleID}-----"

  if (params.gen_org == "human") {
  '''

    !{params.cosmic_annot} \
    -i1 !{params.cosmic} \
    -i2 !{snps_filt} > !{sampleID}_flt_snp_cosmic_annotation.vcf


    java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /snpEff_v4_3/snpEff/snpEff.jar -v -lof !{params.gen_ver} -dataDir !{params.hgvs_data} -noStats !{sampleID}_flt_snp_cosmic_annotation.vcf  > !{sampleID}_SnpEff_snp.vcf


    java -jar /snpEff_v4_3/snpEff/SnpSift.jar dbnsfp -v -db !{params.dbNSFP} -noDownload -a -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF !{sampleID}_SnpEff_snp.vcf > !{sampleID}_full_anno_snp.vcf


    cat !{sampleID}_full_anno_snp.vcf | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_full_anno_snp_onePerLine.vcf

    java -jar /snpEff_v4_3/snpEff/SnpSift.jar extractFields !{sampleID}_full_anno_snp_onePerLine.vcf CHROM POS ID REF ALT QUAL FILTER "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" > !{sampleID}_full_anno_snp.txt


    !{params.cosmic_annot} \
    -i1 !{params.cosmic} \
    -i2 !{indels_filt} > !{sampleID}_flt_indel_cosmic_annotation.vcf


    java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /snpEff_v4_3/snpEff/snpEff.jar -v -lof !{params.gen_ver} -dataDir !{params.hgvs_data} -noStats !{sampleID}_flt_indel_cosmic_annotation.vcf  > !{sampleID}_SnpEff_indel.vcf


    java -jar /snpEff_v4_3/snpEff/SnpSift.jar dbnsfp -v -db !{params.dbNSFP} -noDownload -a -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF !{sampleID}_SnpEff_indel.vcf > !{sampleID}_full_anno_indel.vcf


    cat !{sampleID}_full_anno_indel.vcf | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_full_anno_indel_onePerLine.vcf

    java -jar /snpEff_v4_3/snpEff/SnpSift.jar extractFields !{sampleID}_full_anno_indel_onePerLine.vcf CHROM POS ID REF ALT QUAL FILTER "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" > !{sampleID}_full_anno_indel.txt


    rm -rf !{sampleID}_only_snps.vcf !{sampleID}_only_snps.vcf.idx !{sampleID}_only_snps_filtered.vcf !{sampleID}_only_snps_filtered.vcf.idx !{sampleID}_only_indels.vcf !{sampleID}_only_indels.vcf.idx !{sampleID}_only_indels_filtered.vcf !{sampleID}_only_indels_filtered.vcf.idx

    '''
    }

  else if (params.gen_org == "mouse") {

    '''

    cat !{snps_filt} > !{sampleID}_full_anno_snp.vcf

    cat !{indels_filt} > !{sampleID}_full_anno_indel.vcf

    '''
    }
}

// Step 9: Post Variant Calling Processing, part 2
process post_var_proc2 {
  tag "sampleID"
  label 'short_med_mem'
  label 'gatk3'

  publishDir "${sample_tmpdir}_tmp", pattern: "*snp*", mode: 'copy'

  input:
  tuple sampleID, file(full_annot_snp), file(full_annot_indel) from annot_snp.join(annot_indel)

  output:
  tuple sampleID, file("*snp_indel*.vcf") into snp_indel_combo, dummy_snp_indel_combo
  file("*snp*")

  script:
  log.info "-----Post variant calling processing, part 2 running on ${sampleID}-----"

  if (params.gen_org == "human") {

    """

    java -Djava.io.tmpdir=$TMPDIR -Xmx2g -jar /usr/GenomeAnalysisTK.jar \
    -T CombineVariants \
    --genotypemergeoption UNIQUIFY \
    -R ${params.ref_fa} \
    --variant:SNP ${full_annot_snp} \
    --variant:INDEL ${full_annot_indel} \
    -o ${sampleID}_snp_indel_combined.vcf

    """
    }

  else if (params.gen_org == "mouse") {

    """
    java -Xmx8g -jar /snpEff/snpEff.jar GRCm38.75 \
    -c ${params.snpEff_config} \
    -o gatk \
    -s ${sampleID}_snpeff_summary.html ${full_annot_snp} > ${sampleID}_variants_filtered_dbSNP_snpEff.vcf


    java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
    -R ${params.ref_fa} \
    -T VariantAnnotator \
    -A SnpEff \
    --variant ${full_annot_snp} \
    --snpEffFile ${sampleID}_variants_filtered_dbSNP_snpEff.vcf \
    -L ${full_annot_snp} \
    -o ${sampleID}_snp_indel_variants_filtered_highestsnpEff.vcf


    java -jar /snpEff/SnpSift.jar extractFields ${sampleID}_snp_indel_variants_filtered_highestsnpEff.vcf  CHROM POS REF ALT ID FILTER QUAL FILTER AF SNPEFF_FUNCTIONAL_CLASS SNPEFF_GENE_NAME SNPEFF_AMINO_ACID_CHANGE SNPEFF_EFFECT SNPEFF_TRANSCRIPT_ID > ${sampleID}_snp_indel_variants_filtered_highestsnpEff.txt

    """
    }
}


// Step 10: Aggregate Stats
if (params.gen_org == "human") {
  process agg_stats_hsa {
    tag "sampleID"
    label 'vshort_mem'
    label 'python'

    publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

    input:
    tuple sampleID, file(filter_stats), file(PicardMet), file(CoverageMet) from fq_stats.join(picard_metrics).join(covmet)

    output:
    tuple sampleID, file("*summary_stats.txt")
    tuple sampleID, file("*stats.txt") into dummy_stats_file

    script:
    log.info "-----Generating summary stats file for ${sampleID}-----"

    """

    python ${params.stats_agg} ${sampleID}_summary_stats.txt ${filter_stats} ${PicardMet} ${CoverageMet}


    """
    }
  }

else if (params.gen_org == "mouse") {

  process agg_stats_mmu {
    tag "sampleID"
    label 'vshort_mem'
    label 'python'

    publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

    input:
    tuple sampleID, file(filter_stats), file(AlignMet) from fq_stats.join(covmet)

    output:
    tuple sampleID, file("*summary_stats.txt")
    tuple sampleID, file("*stats.txt") into dummy_stats_file

    script:
    log.info "-----Generating summary stats file for ${sampleID}-----"

    """

    python ${params.stats_agg} ${sampleID}_summary_stats.txt ${filter_stats} ${AlignMet}

    """
    }
  }


// Step 11: Transfer files to output directory
process transfer_files {
  tag "sampleID"
  label 'vshort_mem'

  input:
  tuple sampleID, file(fq_stat) from dummy_fq_stats
  tuple sampleID, file(bwa_mem) from dummy_bwa_mem_out
  tuple sampleID, file(picmet) from dummy_picard_metrics
  tuple sampleID, file(coverage_met) from dummy_covmet
  tuple sampleID, file(raw_var) from dummy_raw_variants
  tuple sampleID, file(gvcf) from dummy_gvcf
  tuple sampleID, file(filtered_var) from dummy_filt_var
  tuple sampleID, file(anno_snp) from dummy_annot_snp
  tuple sampleID, file(anno_indel) from dummy_annot_indel
  tuple sampleID, file(comb_snp_indel) from dummy_snp_indel_combo
  tuple sampleID, file(stat_file) from dummy_stats_file

  script:
  log.info "-----Moving files to output directory for ${sampleID}-----"

  """

  echo "${sample_tmpdir}" > sampletmpdir.txt

  fullpath=`cat sampletmpdir.txt`


  finfullpath=\$(basename \$fullpath)


  mv "\${fullpath}_tmp" "${params.outdir}/\${finfullpath}/"


  """
}

// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Summary ~ ~ ~ ~ ~
workflow.onComplete {
  wfEnd = [:]
  wfEnd['Completed at'] = ${workflow.complete}
  wfEnd['Duration']     = ${workflow.duration}
  wfEnd['Exit status']  = ${workflow.exitStatus}
  wfEnd['Success']      = ${workflow.success}
    if(!workflow.success){
      wfEnd['!!Execution failed'] = ''
      wfEnd['.    Error']   = ${workflow.errorMessage}
      wfEnd['.    Report']  = ${workflow.errorReport}
    }
  Summary.show(wfEnd)
}
