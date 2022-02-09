ef help() {
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
        --_fqPath                 directory containing input fastqs
        --fastqInputs             absolute path to input fastq(s).
        --outdir                  directory where final output will be stored.
        --tmpdir                  directory where temporary files will be held.
        --ref_fa_bwa              absolute path to bwa index .fa file
        --ref_fa                  absolute path to reference genome to be used for alignment
        --dbsnp                   absolute path to dbSNP vcf file
        --filter_trim             script for trimming fastq files
        --read_grp                script to get read group from fastq file
        --gen_org                 organism genome (options are human or mouse)
        --gen_ver                 genome version (options are hg19, hg38, or mm10)
        --min_pct_hq_reads        minimum percentage of HQ reads that remain after trimming
        --stats_agg               script to aggregate stats from multiple steps
        --call_val                threshold to consider site confidently called
        --mismatch_penalty        mismatch penalty for BWA MEM

    Optional:
        --name                  Name for the pipeline run. If not specified Nextflow will
                                automatically generate a random mnemonic.
        --dbNSFP                absolute path to dbNSFP file
        --snpEff_data           absolute path to snpEff data
        --hgvs_data             absolute path to hgvs data; not applicable for mouse sequencing
        --gold_std_indels       absolute path to gold standard indels
        --phase1_1000G          absolute path to 1000 Genomes Phase I indel calls
        --cosmic                absolute path to Cosmic coding and non-coding variants vcf file; not applicable for mouse sequencing
        --cosmic_annot          script to add Cosmic annotations; not applicable for mouse sequencing
        --snpEff_config         absolute path to snpEff config file; not applicable for human sequencing
        --ploidy_val            ploidy value for HaplotypeCaller; not applicable for mouse sequencing
        --V20                   chromosome 20 vcf file; not applicable for mouse sequencing
        --V21                   chromosome 21 vcf file; not applicable for mouse sequencing
        --V22                   chromosome 22 vcf file; not applicable for mouse sequencing



    """.stripIndent()
}
