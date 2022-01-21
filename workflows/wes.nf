#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {BWA_MEM} from '../modules/bwa'

// import groovy files
import Logos
import Helpers

// print logo and set datestamp
logo = new Logo()
println logo.show()
def timestamp = new Date().format("yyMMdd-HHmmss")

// prepare reads channel
if (params.read_type == 'PE'){
  read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}_*${params.extension}",checkExists:true )
}
else if (params.read_type == 'SE'){
  read_ch = Channel.fromFilePairs("${params.fq_path}/*${params.extension}",checkExists:true, size:1 )
}

// main workflow
workflow WES{
  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)
  // Step 2: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)
  // Step 3: BWA-MEM Alignment
  BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq, READ_GROUPS.out.read_groups)
  // Step 4: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.bwa_mem)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.picard_sortsam)
  // Step 5: Variant Pre-Processing - Part 2

}






// Step 1: Qual_Stat
// Step 2a: Get Read Group Information
// Step 2 (3): BWA-MEM Alignment
// Step 4: Variant Preprocessing 1
// Step 5: Variant Pre-Processing 2
process variant_preproc_2 {
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

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

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'picard-1.95.sif'

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

  cpus = 8
  memory = 15.GB
  time = '23:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

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

  cpus = 8
  memory = 15.GB
  time = '23:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

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

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'


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

  cpus = 8
  memory = 10.GB
  time = '08:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

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

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'

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

    cpus = 1
    time = '00:30:00'
    clusterOptions = '-q batch'

    container 'python_2.7.3.sif'

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

    cpus = 1
    time = '00:30:00'
    clusterOptions = '-q batch'

    container 'python_2.7.3.sif'

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

  cpus = 1
  time = '00:30:00'
  clusterOptions = '-q batch'

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
