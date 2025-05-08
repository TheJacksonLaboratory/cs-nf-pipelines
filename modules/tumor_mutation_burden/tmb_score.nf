process TMB_SCORE {
    tag "$sampleID"

    cpus 1
    memory 25.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/snpeff_5.1_r_bedtools_tmb:v3'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*TMB_Score.txt", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), val(tumor_name)
    path(hexcoverage)

    output:
    tuple val(sampleID), file("*TMB_Score.txt"), emit: score
    tuple val(sampleID), file("*HM.tab"), emit: tab
 
    shell:
    '''
    python3 !{projectDir}/bin/wes/allele_depth_min_and_AF_from_ADs.py !{vcf} !{vcf.baseName}.pyfiltered.vcf 15 !{tumor_name}

    java -jar /opt/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' -f !{vcf.baseName}.pyfiltered.vcf | \
    java -jar /opt/snpEff/SnpSift.jar filter --addFilter "strandBias" --rmFilter "PASS" 'FS > 60' | \
    java -jar /opt/snpEff/SnpSift.jar filter --addFilter "lowMQ" --rmFilter "PASS" 'MQ < 40' | \
    java -jar /opt/snpEff/SnpSift.jar filter --addFilter "lowMQRankSum" --rmFilter "PASS" 'MQRankSum < -12.5' | \
    java -jar /opt/snpEff/SnpSift.jar filter --addFilter "lowReadPosRankSum" --rmFilter "PASS" 'ReadPosRankSum < -8' > !{sampleID}_nofalsePositives.tmp2.vcf

    java -jar /opt/snpEff/snpEff.jar eff -v -c !{params.snpEff_config} -lof -canon -hgvs hg38 -noStats !{sampleID}_nofalsePositives.tmp2.vcf > !{sampleID}_all_genes_variants_snpEff.vcf

    java -jar /opt/snpEff/SnpSift.jar dbnsfp -v -db !{params.dbNSFP} -noDownload -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AC,ExAC_AF !{sampleID}_all_genes_variants_snpEff.vcf > !{sampleID}_all_genes_variants_snpeff_snpsift.vcf

    java -jar /opt/snpEff/SnpSift.jar annotate -id !{params.cosmic} !{sampleID}_all_genes_variants_snpeff_snpsift.vcf > !{sampleID}_all_genes_variants_cosmicannotation.vcf

    java -jar /opt/snpEff/SnpSift.jar filter --addFilter "PutativeGermline" --rmFilter "PASS" '(ALT_AF[ANY] >= 90 & dbNSFP_1000Gp3_AF[ANY] >= 0.0095) | ((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_1000Gp3_AF[ANY] >= 0.0095)| (ALT_AF[ANY] >= 90 & dbNSFP_ESP6500_AA_AF[ANY] >= 0.0095)| ((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ESP6500_AA_AF[ANY] >= 0.0095)|(ALT_AF[ANY] >= 90 & dbNSFP_ESP6500_EA_AF[ANY] >= 0.0095)|((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ESP6500_EA_AF[ANY] >= 0.0095)|(ALT_AF[ANY] >= 90 & dbNSFP_ExAC_AF[ANY] >= 0.0095)|((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ExAC_AF[ANY] >= 0.0095)' !{sampleID}_all_genes_variants_cosmicannotation.vcf > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag.vcf

    cat !{sampleID}_all_genes_variants_cosmicannotation_germlineflag.vcf | /opt/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline.vcf

    java -jar /opt/snpEff/SnpSift.jar extractFields -noDownload -noLog !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline.vcf "CHROM" "POS" "REF" "ALT" "ID" "FILTER" "DP" "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_ExAC_AC" "dbNSFP_ExAC_AF" > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline_final.tab

    cat !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline_final.tab | grep "HIGH\\|MODERATE" | awk -F '\\t' '{ if($5 !~ "rs" && ($6==""||$6=="PASS"||$6==".")) print $1,$2-1,$2 }'| sort| uniq| tr ' ' '\\t' > count2; echo "chr\tstart\tend\tLength\tHM" > !{sampleID}_somatic.HM.tab; bedtools coverage -a !{hexcoverage} -b count2 | cut -f 1-5 >> !{sampleID}_somatic.HM.tab

    Rscript !{projectDir}/bin/wes/TMB_calc.R !{sampleID}_somatic.HM.tab  !{sampleID}_somatic_TMB_Score.txt

    '''
}
