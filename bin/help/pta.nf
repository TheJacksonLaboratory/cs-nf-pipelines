def help(){
if (params.gen_org=='human')
println '''
Parameter | Default | Description

The following are human specific parameters. To see help for mouse, add `--gen_org mouse` to your command. 

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--csv_input | /<FILE_PATH> | CSV delimited sample sheet that controls how samples are processed. The required input header is: patient,sex,status,sampleID,lane,fastq_1,fastq_2. See the repository wiki (https://github.com/TheJacksonLaboratory/cs-nf-pipelines/wiki) for additional information. 

--xenome_prefix | /projects/compsci/omics_share/human/GRCh38/supporting_files/xenome/trans_human_GRCh38_84_NOD_based_on_mm10_k25| Xenome index for deconvolution of human and mouse reads. Used when `--pdx` is run. 
--pdx | false | Options: false, true. If specified, 'Xenome' is run on reads to deconvolute human and mouse reads. Human only reads are used in analysis. 

--deduplicate_reads | false | Options: false, true. If specified, run bbmap clumpify on input reads. Clumpify will deduplicate reads prior to trimming. This can help with mapping and downstream steps when analyzing high coverage WGS data. 

--coverage_cap | null | If an integer value is specified, jvarkit 'Biostar154220' is used to cap coverage at the that value. See: http://lindenb.github.io/jvarkit/Biostar154220.html
--primary_chrom_bed | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/intervals/Homo_sapiens_assembly38.primary_chrom.bed' | A bed file containing the primary chromsomes with positions. Used in limiting jvarkit 'Biostar154220' to those regions with expected coverage.

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.

--ref_fa | '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta' | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. JAX users should not change this parameter.
--ref_fa_indices | '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta' | Pre-compiled BWA index files. JAX users should not change this parameter.

--ref_fa_dict | '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.dict' | FASTA dictonary file. JAX users should not change this parameter. 
--combined_reference_set | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/combined_ref_set/Homo_sapiens_assembly38.fasta' | Several tools (GRIDSS, SVABA) requires reference and bwa index files in same directory. Links used within this directory to avoid duplication of fasta and bwa indicies. See note in directory. 

--mismatch_penalty | -B 8 | The BWA penalty for a mismatch.

--gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Used in GATK BaseRecalibrator and variant tranche recalibration derived from the GATK resource bundle. JAX users should not change this parameter.
--phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Used in GATK BaseRecalibrator derived from the GATK resource bundle. JAX users should not change this parameter.
--dbSNP | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz' | Used in variant annotation, GATK BaseRecalibrator, variant tranche recalibration, and by SVABA. JAX users should not change this parameter.
--dbSNP_index | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz.tbi' | Index associated with the dbsnp file. 

--chrom_contigs | '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.primaryChr.contig_list' | Contig list used for scatter / gather in calling and annotation. 
--chrom_intervals | '/projects/omics_share/human/GRCh38/genome/annotation/intervals/hg38_calling_intervals/' | Chromosome intervals used for scatter gather in calling. 

--call_val | 50 | The minimum phred-scaled confidence threshold at which variants should be called.
--ploidy_val | '-ploidy 2' | Sample ploidy used by Haplotypecaller in germline small variant calling. 

--excludeIntervalList | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/intervals/hg38_haplotypeCaller_skip.interval_list' | Germline caller exclusion list. 
--hapmap | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/snps_indels/hapmap_3.3.hg38.vcf.gz' | variant tranche recalibration requirement derived from the GATK resource bundle. 
--omni | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_omni2.5.hg38.vcf.gz' | variant tranche recalibration requirement derived from GATK resource bundle. 

--pon_bed | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/filtering/WGS_1000g_GRCh38.pon.bed' | Panel of normal samples used in in snp and indel filtering.
--intervalListBed | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/filtering/SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla.interval_list.bed' | This file is used to extract small variants in non-exonic regions. Such calls are then attempted to be recovered via Lancet calls.

--lancet_beds_directory | '/projects/omics_share/human/GRCh38/genome/annotation/intervals/lancet_chr_beds/' | Lancet interval bed files used in calling by that tool. 

--mappability_directory | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/intervals/mappability' | Bicseq2 input requirement. Derived from the tool developer resource pack. 
--bicseq2_chromList | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/configs/sampleId.bicseq2.config' | Bicseq2 config requirement. Derived from the tool developer resource pack.
--bicseq2_no_scaling | false | false: estimate 'lamda' smoothing factor from data for CNV profile calling. true: Use standard 'lamda | 4' smoothing for CNV profile calling. If BicSeq2 fails with an error, set this parameter to 'true'.

--gripss_pon | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/gripss_pon' | Panel of normal files for Gripss SV call filering. Provided by the tool developer resource pack.

--callRegions | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/intervals/GRCh38.callregions.bed.gz' | Manta calling regions. Provided by the tool developer resource pack.

--strelka_config | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/configs/configureStrelkaSomaticWorkflow.py.ini' | Strelka input configuration. Provided by the tool developer resource pack.

--msisensor_model | '/projects/compsci/omics_share/human/GRCh38/supporting_files/msisensor2/models_hg38' | Model files for MSI calling via MSIsensor2. Provided by the tool developer resource pack.

--vep_cache_directory | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/vep_data' | VEP annotation cache. Cache provided is for Ensembl v109.  
--vep_fasta | '/projects/compsci/omics_share/human/GRCh38/genome/sequence/ensembl/GRCh38.p13/Homo_sapiens.GRCh38.dna.primary_assembly.fa' | VEP requires an ensembl based fasta. GRCh38.p13 is used for v97-v109.  

--cosmic_cgc | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/function/cancer_gene_census_v97.csv' | COSMIC Cancer Gene Census annotation file. Index for file required within same location. 
--cosmic_cancer_resistance_muts | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/function/CosmicResistanceMutations.tsv.gz' | COSMIC Resistance Mutations file. Index for file required within same location. 

--ensembl_entrez | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/GRCh39.p13_ensemblv109_entrez_id_map.csv' | Ensembl to Entrez gene ID to HGNC symbol mapping file. used in somatic vcf finalization.

--germline_filtering_vcf | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/filtering/gnomad-and-ALL_GRCh38_sites.20170504.normalized.modified.PASS.vcf.gz' | Germline reference file used in SNV call filtering. 
--cytoband | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/GRCh38.cytoBand.UCSC.chr.sorted.txt' | File used in bicseq2 annotations
--dgv | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/DGV.GRCh38_hg38_variants_2020-02-25.bed' | File used in bicseq2 annotations
--thousandG | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/1KGP.CNV.GRCh38.canvas.merged.bed' | File used in bicseq2 annotations
--cosmicUniqueBed | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/CosmicCompleteCNA_uniqIntervals.bed' | File used in bicseq2 annotations
--cancerCensusBed | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/cancer_gene_census.GRCh38-v92.bed' | File used in bicseq2 annotations and SV annotation. 
--ensemblUniqueBed | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/ensembl_genes_unique_sorted.final.v93.chr.sorted.bed' | File used in bicseq2 annotations and SV annotation.   
--gap | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/GRCh38.gap.UCSC.annotated.chr.sorted.bed' | File used in SV annotation.
--dgvBedpe | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/DGV.GRCh38_hg38_variants_2020-02-25.bedpe' | File used in SV annotation.
--thousandGVcf | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/1KGP.pruned_wAFs.PASS_and_MULTIALLELIC_Mosaic.GRCh38.vcf' | File used in SV annotation.
--svPon | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/1000G-SV-PON.survivor-merged.GRCh38.filtered.bedpe' | File used in SV annotation.
--cosmicBedPe | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/cosmic-sv-GRCh38-v92.bedpe' | File used in SV annotation.

--na12878_bam | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/NA12878/NA12878_realigned_BQSR.bam' | NA12878 BAM file. Used in un-paired sample analysis. 
--na12878_bai | '/projects/compsci/omics_share/human/GRCh38/supporting_files/PTA_inputs/NA12878/NA12878_realigned_BQSR.bai' | NA12878 BAM index file. Used in un-paired sample analysis. 
--na12878_sampleName | 'ERR194147_1.fastq.gz_filtered_trimmed' | NA12878 sample name within the NA12878 BAM file. 

--read_type | PE | Only 'PE' is accepted for this workflow. 

'''
else
println '''

The following are mouse specific parameters. To see help for mouse, add `--gen_org human` to your command. 

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--csv_input | /<FILE_PATH> | CSV delimited sample sheet that controls how samples are processed. The required input header is: patient,sex,status,sampleID,lane,fastq_1,fastq_2. See the repository wiki (https://github.com/TheJacksonLaboratory/cs-nf-pipelines/wiki) for additional information. 

--deduplicate_reads | false | Options: false, true. If specified, run bbmap clumpify on input reads. Clumpify will deduplicate reads prior to trimming. This can help with mapping and downstream steps when analyzing high coverage WGS data. 

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.

--ref_fa | '/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/ensembl/GRCm39.p0/Mus_musculus.GRCm39.dna.primary_assembly.fa' | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. JAX users should not change this parameter.
--ref_fa_indices | '/projects/compsci/omics_share/mouse/GRCm39/genome/indices/ensembl/v105/bwa/Mus_musculus.GRCm39.dna.primary_assembly.fa' | Pre-compiled BWA index files. JAX users should not change this parameter.
--ref_fa_dict | '/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/ensembl/GRCm39.p0/Mus_musculus.GRCm39.dna.primary_assembly.dict' | FASTA dictonary file. JAX users should not change this parameter. 
--combined_reference_set | '/projects/compsci/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/combined_ref_set/Mus_musculus.GRCm39.dna.primary_assembly.fa' | Several tools (GRIDSS, SVABA) requires reference and bwa index files in same directory. Links used within this directory to avoid duplication of fasta and bwa indicies. See note in directory. 

--mismatch_penalty | -B 8 | The BWA penalty for a mismatch.

--dbSNP | '/projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/GCA_000001635.9_current_ids.vcf.gz' | Used in variant annotation and by SVABA. JAX users should not change this parameter.
--dbSNP_index | '/projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/GCA_000001635.9_current_ids.vcf.gz.tbi' | Index associated with the dbsnp file. 

--chrom_contigs | '/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/ensembl/GRCm39.p0/Mus_musculus.GRCm39.dna.primary_assembly.primaryChr.contig_list' | Contig list used for scatter / gather in calling and annotation. 
--chrom_intervals | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/GRCm39_calling_intervals/' | Chromosome intervals used for scatter gather in calling. 

--call_val | 50 | The minimum phred-scaled confidence threshold at which variants should be called.
--ploidy_val | '-ploidy 2' | Sample ploidy used by Haplotypecaller in germline small variant calling. 

--excludeIntervalList | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/mm39.excluderanges.interval_list' | Germline caller exclusion list. 

--intervalListBed | '/projects/compsci/omics_share/mouse/GRCm39/supporting_files/capture_kit_files/agilent/v2/S32371113_mouse_exon_V2.mm39.bare.bed' | This file is used to extract small variants in non-exonic regions. Such calls are then attempted to be recovered via Lancet calls.

--lancet_beds_directory | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/lancet_chr_beds/' | Lancet interval bed files used in calling by that tool. 

--delly_exclusion | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/GRCm39_gap_delly_exclusion.txt' | Delly CNV calling exclusion list. 
--delly_mappability | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/mappability/GRCm39.p0.map.gz' | Delly CNV calling mappability file. 
--cnv_window | 10000 | Delly CNV calling read depth window size. Default value is tool default. This parameter is included for testing purposes only.
--cnv_min_size | 10000 | Delly CNV classification minimum CNV size. Default value is tool default. This parameter is included for testing purposes only.
--cnv_germline_prob | 0.00100000005 | Delly CNV classification germline probability. Default value is tool default. This parameter is included for testing purposes only. 

--callRegions | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/GRCm39.callregions.bed.gz' | Manta calling regions. Provided by the tool developer resource pack.

--strelka_config | '/projects/compsci/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/configs/configureStrelkaSomaticWorkflow.py.ini' | Strelka input configuration. Provided by the tool developer resource pack.

--vep_cache_directory | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/vep_data' | VEP annotation cache. Cache provided is for Ensembl v109.  
--vep_fasta | '/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/ensembl/GRCm39.p0/Mus_musculus.GRCm39.dna.primary_assembly.fa'  | VEP requires an ensembl based fasta. GRCh38.p13 is used for v97-v109.  

--cytoband | '/projects/compsci/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/annotations/GRCm38.liftedTo.GRCm39.cytoBand.UCSC.chr.sorted.bed' | Cytoband file used in CNV annotations. Derived from UCSC table, lifted from GRCm38 to GRCm39. 
--known_del | '/projects/omics_share/mouse/GRCm39/genome/annotation/struct_vars/ferraj_2023_inv_ins_del/variants_freeze5_sv_sym_DEL_mm39_sorted.bed' | Used in SV annotation, and filtering.
--known_ins | '/projects/omics_share/mouse/GRCm39/genome/annotation/struct_vars/ferraj_2023_inv_ins_del/variants_freeze5_sv_sym_INS_mm39_sorted.bed' | Used in SV annotation, and filtering.
--known_inv | '/projects/omics_share/mouse/GRCm39/genome/annotation/struct_vars/ferraj_2023_inv_ins_del/variants_freeze5_sv_sym_INV_mm39_sorted.bed' | Used in SV annotation, and filtering.

--ensemblUniqueBed | '/projects/compsci/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/annotations/ensembl_genes_unique_sorted.final.v110.chr.sorted.bed' | Ensembl gene IDs and positions, used in SV and CNV annotations. 

--gap | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/GRCm39_gap.bed' | File used in SV annotation. From UCSC table browser.
--exclude_list | '/projects/compsci/omics_share/mouse/GRCm39/genome/annotation/intervals/mm39.excluderanges_cleaned.bed' | File used in SV annotation. From: https://dozmorovlab.github.io/excluderanges/.

--proxy_normal_bam | '/projects/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/C57L_J/C57L_J_dedup.bam' | C57L_J BAM file. Used in un-paired sample analysis. 
--proxy_normal_bai | '/projects/omics_share/mouse/GRCm39/supporting_files/PTA_inputs/C57L_J/C57L_J_dedup.bam.bai' | C57L_J BAM index file. Used in un-paired sample analysis. 
--proxy_normal_sampleName | 'C57L_J' | C57L_J sample name within the C57L_J BAM file. 

--read_type | PE | Only 'PE' is accepted for this workflow. 
'''
}
