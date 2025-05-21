process VEP_SOMATIC {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'ensemblorg/ensembl-vep:release_109.3'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(vcf), file(idx), val(meta), val(normal_name), val(tumor_name)

    output:
    tuple val(sampleID), file("*_vep_annotated.vcf"), val(meta), val(normal_name), val(tumor_name), emit: vcf

    script:

    """
    vep \
    --input_file ${vcf} \
    --output_file ${sampleID}_somatic_vep_annotated.vcf \
    --fork ${task.cpus} \
    --buffer_size 50000 \
    --format vcf \
    --no_stats \
    --no_escape \
    --offline \
    --assembly GRCh38 \
    --cache \
    --dir_cache ${params.vep_cache_directory} \
    --refseq \
    --max_af \
    --af \
    --af_1kg \
    --af_gnomad \
    --exclude_predicted \
    --fasta ${params.vep_fasta} \
    --symbol \
    --hgvs \
    --check_existing \
    --vcf \
    --pick_allele_gene \
    --dir_plugins ${params.vep_cache_directory}/Plugins \
    --plugin dbscSNV,${params.vep_cache_directory}/Plugins/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz \
    --plugin MaxEntScan,${params.vep_cache_directory}/Plugins/maxentscan \
    --plugin dbNSFP,${params.vep_cache_directory}/Plugins/dbNSFP/dbNSFP4.3a_grch38.gz,${params.vep_cache_directory}/Plugins/dbNSFP_replacement_logic,REVEL_score,SIFT_pred,SIFT4G_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaSVM_pred,PrimateAI_pred,fathmm-MKL_coding_pred,GERP++_RS,phyloP100way_vertebrate,CADD_phred,Polyphen2_HVAR_pred \
    --custom ${params.vep_cache_directory}/annotations/COSMIC_v97/CosmicCodingMuts.vcf.gz,CosmicCoding,vcf,exact,0,GENOMIC_ID,LEGACY_ID,CNT,CDS,AA \
    --custom ${params.vep_cache_directory}/annotations/COSMIC_v97/CosmicNonCodingVariants.normal.vcf.gz,CosmicNonCoding,vcf,exact,0,GENOMIC_ID,LEGACY_ID,CNT,CDS,AA \
    --custom ${params.vep_cache_directory}/annotations/04142020_NYGC_samples.vcf.gz,NYGC,vcf,exact,0,AF,Samples,AC_Het,AC_Hom \
    --custom ${params.vep_cache_directory}/annotations/clinvar.vep.vcf.gz,CLN_Overlap,vcf,overlap,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
    --custom ${params.vep_cache_directory}/annotations/clinvar.vep.vcf.gz,CLN_Exact,vcf,exact,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
    --custom ${params.vep_cache_directory}/annotations/gnomad_exomes_subset_final.vcf.gz,GnomadExomes,vcf,exact,0,AF,nhomalt \
    --custom ${params.vep_cache_directory}/annotations/gnomad_genomes_subset_final.vcf.gz,GnomadGenomes,vcf,exact,0,AF,nhomalt \
    --custom ${params.vep_cache_directory}/annotations/chd_genes.vcf.gz,CHD_GENES,vcf,overlap,0,GENE \
    --custom ${params.vep_cache_directory}/annotations/chd_evolving.vcf.gz,CHD_EVOLVING,vcf,overlap,0,GENE \
    --custom ${params.vep_cache_directory}/annotations/chd_whitelist.vcf.gz,chd_whitelist,vcf,overlap,0,END \
    --custom ${params.vep_cache_directory}/annotations/deep_intronic_whitelist_08132020.vcf.gz,INTRONIC,vcf,exact,0,INTRONIC \
    --custom ${params.vep_cache_directory}/annotations/clinvar_deep_intronics_09012020.vcf.gz,CLINVAR_INTRONIC,vcf,exact,0,INTRONIC \
    --custom ${params.vep_cache_directory}/annotations/mastermind_cited_variants_reference-2021.01.02-grch38_fixed-contigs.vcf.gz,mm,vcf,exact,0,GENE,HGVSG,MMCNT1,MMCNT2,MMCNT3,MMID3,MMURI3 \
    --custom ${params.vep_cache_directory}/annotations/spliceai_scores.hg38.sorted.vcf.gz,SPLICEAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL \
    --custom ${params.vep_cache_directory}/annotations/pli_hg38.vcf.gz,PLI,vcf,overlap,0,pLI,mis_z \
    --custom ${params.vep_cache_directory}/annotations/domino_genes_38.vcf.gz,Domino,vcf,overlap,0,Domino_Score \
    --custom ${params.vep_cache_directory}/annotations/ar_extended.vcf.gz,AR,vcf,overlap,0,AR_GENE \
    --custom ${params.vep_cache_directory}/annotations/ACMG59_2017-09-28.vcf.gz,ACMG59,vcf,overlap,0,GENE,DISEASE \
    --custom ${params.vep_cache_directory}/annotations/dials_genes_b38.vcf.gz,DIALS,vcf,overlap,0,DIALS_GENE \
    --custom ${params.vep_cache_directory}/annotations/pgx_vep.vcf.gz,PGx,vcf,exact,0,pgx_rsid \
    --custom ${params.vep_cache_directory}/annotations/sema4_immuno_genes_b38.vcf.gz,IMMUNO,vcf,overlap,0,IMMUNO_Gene \
    --custom ${params.vep_cache_directory}/annotations/sema4_neuro_genes_b38.vcf.gz,NEURO,vcf,overlap,0,NEURO_Gene \
    --custom ${params.vep_cache_directory}/annotations/sema4_cardio_genes_b38.vcf.gz,CARDIO,vcf,overlap,0,CARDIO_Gene \
    --custom ${params.vep_cache_directory}/annotations/nygc_curation_b38.vcf.gz,N19,vcf,overlap,0,NYGC_CUR \
    --custom ${params.vep_cache_directory}/annotations/nygc_reported_variants_b38.vcf.gz,R19,vcf,overlap,0,NYGC_REPORTED_SAMPLE,NYGC_CLASS,NYGC_DISEASE
    """
}

// NOTE: Many of the resources are hard coded based on those provided in: 
//       https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/annotate/variantEffectPredictor.wdl
//       https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/config/fasta_references.json
//       For VEP cache, dbNSFP and dbscSNV resources were rebuild using VEPv108 and as noted below. 

// VEP Cache setup: 

// singularity pull --name vep.sif docker://ensemblorg/ensembl-vep:release_109.3
// singularity exec vep.sif INSTALL.pl -c /PATH_TO_VEP/vep -a cfp -s homo_sapiens_refseq -y GRCh38 -g dbNSFP,dbscSNV,MaxEntScan

// In the plugin directory: 

// dbNSFP: 
// wget https://dbnsfp.s3.amazonaws.com/dbNSFP4.3a.zip
// unzip dbNSFP4.3a.zip
// zcat dbNSFP4.3a_variant.chr1.gz | head -n1 > h
// mkdir temp
// zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | sort -T temp -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.3a.gz
// tabix -s 1 -b 2 -e 2 dbNSFP4.3a.gz
// rm -rf temp dbNSFP4.3a_variant.chr* h dbNSFP4.3_gene.gz dbNSFP4.3_gene.complete.gz dbNSFP4.3a.zip search_dbNSFP43a.readme.pdf search_dbNSFP43a.class search_dbNSFP43a.jar


// dbscSNV:
// wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip
// unzip dbscSNV1.1.zip
// head -n1 dbscSNV1.1.chr1 > h2
// mkdir temp2
// cat dbscSNV1.1.chr* | grep -v ^chr | sort -T temp2 -k5,5 -k6,6n | cat h2 - | awk '$5 != "."' | bgzip -c > dbscSNV1.1_GRCh38.txt.gz
// tabix -s 5 -b 6 -e 6 -c c dbscSNV1.1_GRCh38.txt.gz
// rm dbscSNV1.1.chr* dbscSNV1.1.zip h2

// wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
// gunzip fordownload.tar.gz
// tar -xvf fordownload.tar
// mkdir maxentscan
// mv fordownload ./maxentscan/


// COSMIC: 
// mkdir COSMIC_v97
// echo "<EMAIL>@jax.org:<PASSWORD>" | base64
// <base64 string>
// curl -H "Authorization: Basic ADD AUTHORIZATION" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/VCF/CosmicCodingMuts.vcf.gz
// the above command provides a URL for curl download
// curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v97/VCF/CosmicCodingMuts.vcf.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1672843877&Signature=TSsIiQodqoKS5skE1ziS49zEWSU%3D" --output CosmicCodingMuts.vcf.gz

// curl -H "Authorization: Basic ADD AUTHORIZATION" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/VCF/CosmicNonCodingVariants.normal.vcf.gz
// the above command provides a URL for curl download
// curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v97/VCF/CosmicNonCodingVariants.normal.vcf.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1672844121&Signature=4JkeRizNMg0pv%2FChw4QAl268dVw%3D" --output CosmicNonCodingVariants.normal.vcf.gz


// ALL REMAINING ANNOTATIONS: 
// Note: remaining annotations come from: gs://nygc-resources-public/ensembl_vep/annotations.tar.gz