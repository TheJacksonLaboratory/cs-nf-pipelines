#!/usr.bin.env nextflow
/*
========================================================================
                     mmrSVD pipeline
========================================================================

The Mouse Mutant Resource Structural Variant Detection (mmrSVD) pipeline,
developed by Brian Sanderson (brian.sanderson@jax.org) and 
Mike Lloyd (mike.lloyd@jax.org) for The Jackson Laboratory.

*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow -c /path/to/params.config run /path/to/mmrSVD/main.nf -profile slurm,singularity --genome mm10

    """
}

// Show help message
if (params.help) exit 0, helpMessage()

// fasta can be either given as a genome name in iGenomes or as a fasta file
// params.fasta will be changed only if it's not previously defined
params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : "null"
ch_fastq1 = params.fastq1 ? Channel.value(file(params.fastq1)) : null
ch_fastq2 = params.fastq2 ? Channel.value(file(params.fastq2)) : null
params.surv_dist = 1000
params.surv_supp = 1
params.surv_type = 1
params.surv_strand = 1
params.surv_min = 30
def sample_name = params.names
def abs_outdir = params.outdir
params.threads = 8
params.keep_intermediate = false

// ========================================= Long Read Alignment Processes =====

if (params.seqmode == 'pacbio') {
	
	ch_pbsvTandem = params.pbsv_tandemrepeats? Channel.value(file(params.pbsv_tandemrepeats)) : null

	// CALLERS

	// Sniffels
	
	// https://github.com/fritzsedlazeck/Sniffles
	// Note: Sniffles prefers NGMLR based alignment
	
	// https://github.com/philres/ngmlr
	
	// ALIGN
	// ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam
	
	process NGMLRmap{
		label 'cpus_24'
		label 'ngmlr'
		stageInMode 'copy'
		input:
			val name_string from params.names
			file fasta from ch_fasta
			file fq1 from ch_fastq1
		output:
			file "${name_string}.sam" into ngmlr_sam
		script:
		"""
		ngmlr -t ${task.cpus} --bam-fix -r ${fasta} -q ${fq1} -o ${name_string}.sam
		"""
	}	
	
	
	// SORT BAM
		
	process NGMLRsort{
	  publishDir params.outdir, mode:'copy'
	  label 'cpus_8'
	  label 'samtools_1_9'
	  input:
	  	val name_string from params.names
		file sam from ngmlr_sam
	  output:
		file "${name_string}.ngmlr.aligned.bam" into ch_bam_map
		file "${name_string}.ngmlr.aligned.bam.bai" into ch_bam_index
	  script:
	  """
	  samtools sort --threads ${task.cpus} -m 30${sam} > ${name_string}.ngmlr.aligned.bam
	  samtools index ${name_string}.ngmlr.aligned.bam
	  """
	}
	
	// For Oxford Nanopore run:

	// ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam -x ont
	
	
	// CALL WITH SNIFFLES 
	
	process sniffles {
		publishDir params.outdir, mode:'copy'
		label 'cpus_8'
		label 'sniffles'
		input:
			file bam from ch_bam_map
		output:
			file "*.vcf" into sniffles_vcf
			path(vcf_path) into vcf_sniffles_path
		script:
			  """
			  sniffles -m ${bam} -v sniffles_calls.vcf
			  echo ${params.outdir}/sniffles_calls.vcf > vcf_path # for later merging
			  """
	}
	
	
	// CALL WITH SVIM 
	
	process svim {
		publishDir params.outdir, mode:'copy'
		label 'cpus_8'
		label 'svim'
		input:
			val name_string from params.names
			file bam from ch_bam_map
			file fasta from ch_fasta
		output:
			path "${name_string}_svim" into svim_output
			file "svim_variants.vcf" into svim_vcf
		script:
			  """
			  svim alignment ${name_string}_svim ${bam} ${fasta}
			  cp ${name_string}_svim/variants.vcf svim_variants.vcf
			  """
	}
	
	
	// CALL WITH CUTESV 
	
	if (params.pbmode == 'ccs') {
		process cutesv_css {
			publishDir params.outdir, mode:'copy'
			label 'cpus_8'
			label 'cutesv'
			input:
				file bam from ch_bam_map
				file bai from ch_bam_index
				file fasta from ch_fasta
			output:
				file "cutesv_calls.vcf" into cutesv_css_vcf
			script:
				  """
				  cuteSV ${bam} ${fasta} cutesv_calls.vcf ./ \
				  --threads ${task.cpus} \
                  --max_cluster_bias_INS 1000 \
                  --diff_ratio_merging_INS 0.9 \
                  --max_cluster_bias_DEL 1000 \
                  --diff_ratio_merging_DEL 0.5
				  """
		}
	}
	
	if (params.pbmode == 'clr') {
		process cutesv_clr {
			publishDir params.outdir, mode:'copy'
			label 'cpus_8'
			label 'cutesv'
			input:
				file bam from ch_bam_map
				file bai from ch_bam_index
				file fasta from ch_fasta
			output:
				file "cutesv_calls.vcf" into cutesv_clr_vcf
			script:
				  """
				  cuteSV ${bam} ${fasta} cutesv_calls.vcf ./ \
				  --threads ${task.cpus} \
                  --max_cluster_bias_INS 100 \
                  --diff_ratio_merging_INS 0.3 \
                  --max_cluster_bias_DEL 200 \
                  --diff_ratio_merging_DEL 0.5
				  """
		}
	}

	// PBSV
	
	// NOTE: PBSV prefers PBMM2 based alignments
	
	// https://github.com/PacificBiosciences/pbmm2

	// BUILD INDEX

	process BuildPBMM2index {
		tag "${fasta}"
		label 'cpus_8'
		label 'pbmm2'
		stageInMode 'copy'
		input:
			file fasta from ch_fasta
		output:
			file("${fasta.baseName}.mmi") into pbmm2_built
		script:
		"""
		pbmm2 index ${fasta} ${fasta.baseName}.mmi	
		"""
	}
	
	// ALIGN CCS FASTQ DATA
	if (params.pbmode == 'ccs') {
		process PBMM2fastqMap_css{
			label 'cpus_8'
			label 'pbmm2'
			publishDir params.outdir, mode:'copy'
			input:
				val name_string from params.names
				file mmi from pbmm2_built
				file fq1 from ch_fastq1

			output:
				file "${name_string}.pbmm2.aligned.bam" into pbmm2_bam_css
			script:
			"""
			pbmm2 align ${mmi} ${fq1} ${name_string}.pbmm2.aligned.bam --preset CCS --sort -j ${task.cpus}
			"""
		}	
		// NOTE: removed `--rg \$(cat $rgr)`. It wasn't working. 
		// NOTE: removed `file rgr from readgroup` from inputs
	
		pbmm2_bam = pbmm2_bam_css 
	
	}
	
	if (params.pbmode == 'clr') {
		process PBMM2fastqMap_clr{
			label 'cpus_8'
			label 'pbmm2'
			publishDir params.outdir, mode:'copy'
			input:
				val name_string from params.names
				file mmi from pbmm2_built
				file fq1 from ch_fastq1

			output:
				file "${name_string}.pbmm2.aligned.bam" into pbmm2_bam_clr
			script:
			"""
			pbmm2 align ${mmi} ${fq1} ${name_string}.pbmm2.aligned.bam --median-filter --sort -j ${task.cpus}
			"""
		}	
		
		pbmm2_bam = pbmm2_bam_clr 
		
	}
	
	// CALL WITH PBSV
	
	// Discover signatures of structural variation
	
	
	if (ch_pbsvTandem == null) {
		process pbsv_discovery_tandem{
			label 'cpus_8'
			label 'pbsv'
			input:
				file bam from pbmm2_bam
			output:
				file "*.svsig.gz" into pbsv_svsig_tandem
			when: !(ch_pbsvTandem)
			script:
			"""
			/usr/bin/env bash ${projectDir}/bin/pbsv_tandem.sh ${bam} ${bam.baseName}.svsig.gz
			"""
		}	
	
		pbsv_svsig = pbsv_svsig_tandem
	}
	
	if (ch_pbsvTandem != null) {
	
	
		process pbsv_discovery_no_tandem{
			label 'cpus_8'
			label 'pbsv'
			input:
				file bam from pbmm2_bam
				file tandem from ch_pbsvTandem
			output:
				file "*.svsig.gz" into pbsv_svsig_no_tandem
			when: ch_pbsvTandem
			script:
			"""
			pbsv discover ${bam} ${bam.baseName}.svsig.gz
			"""
		}
		pbsv_svsig = pbsv_svsig_no_tandem
	}


	if (params.pbmode == 'ccs') {
		process pbsv_call_ccs{
			publishDir params.outdir, mode:'copy'
			label 'cpus_8'
			label 'pbsv'
			stageInMode 'copy'
			input:
				file svsig from pbsv_svsig
				file fasta from ch_fasta
			output:
				file "*.vcf" into pbsv_vcf_ccs
				path(vcf_path) into vcf_pbsv_path
			"""
			pbsv call --ccs ${fasta} ${svsig} pbsv_calls.vcf
			echo ${params.outdir}/pbsv_calls.vcf > vcf_path # for later merging
			"""
		}
	}	


	if (params.pbmode == 'clr') {
		process pbsv_call_clr{
			publishDir params.outdir, mode:'copy'
			label 'cpus_8'
			label 'pbsv'
			stageInMode 'copy'
			input:
				file svsig from pbsv_svsig
				file fasta from ch_fasta
			output:
				file "*.vcf" into pbsv_vcf_clr
				path(vcf_path) into vcf_pbsv_path
			"""
			pbsv call ${fasta} ${svsig} pbsv_calls.vcf
			echo ${params.outdir}/pbsv_calls.vcf > vcf_path # for later merging
			"""
		}
	}	
	
	process prep_vcf_list{
    label 'tiny_job'
    input:
        file "pbsv_calls.vcf" from pbsv_vcf_ccs
        file "sniffles_calls.vcf" from sniffles_vcf
    output:
        file "vcf_list" into ch_vcf
    script:
    """
    ls pbsv_calls.vcf >> vcf_list
    ls sniffles_calls.vcf >> vcf_list
    """
	}


	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Collect all VCF file Paths  ~~~~~
	// NOTE: merged extension contains 'BDLM' => Breakdancer + Delly + Lumpy + Manta
	vcf_pbsv_path.concat(
	vcf_sniffles_path
	)
	.collectFile(name: "sample.vcfs.txt", sort: false)
	.set { sample_vcfs_paths }

	/*process survivor{
		label 'long_himem'
		label 'survivor'
		input:
			file "vcf_list" from ch_vcf
			file "pbsv_calls.vcf" from pbsv_vcf_ccs
			file "sniffles_calls.vcf" from sniffles_vcf       
			val name_string from params.names
			val surv_dist from params.surv_dist
			val surv_supp from params.surv_supp
			val surv_type from params.surv_type
			val surv_strand from params.surv_strand
			val surv_min from params.surv_min
		output:
			file "${name_string}.merged.vcf" into ch_merged_vcf
		script:
		"""
		SURVIVOR merge vcf_list ${surv_dist} ${surv_supp} ${surv_type} ${surv_strand} 0 ${surv_min} ${name_string}.merged.vcf
		"""
	} */
	/*
	process annotate_sv{
		label 'mid_job'
		input:
			file merged_vcf from ch_merged_vcf
			val name_string from params.names
		output:
			file "${name_string}.merged.overlap.annotated.txt" into ch_annot
		script:
		"""
		/usr/bin/env bash ${projectDir}/bin/surv_annot.sh ${name_string} ${merged_vcf}
		"""
	}

	process summarize_sv{
		label 'mid_job'
		label 'pyvcf'
		input:
			file merged_vcf from ch_merged_vcf
			val name_string from params.names
		output:
			file "${name_string}.survivor_summary.csv" into ch_summary
		script:
		"""
		/usr/bin/env python ${projectDir}/bin/sv_to_table.py -v ${merged_vcf} -o ${name_string}.survivor_summary.csv
		"""
	}

	process prep_beds{
		label 'short_himem'
		label 'tidyverse'
		input:
			file annot from ch_annot
			file summary from ch_summary
			val name_string from params.names
		output:
			file "${name_string}.ins.bed" into ch_ins
			file "${name_string}.inv.bed" into ch_inv
			file "${name_string}.del.bed" into ch_del
			file "${name_string}.dup.bed" into ch_dup
			file "${name_string}.tra.bed" into ch_tra
		script:
		"""
		/usr/bin/env Rscript ${projectDir}/bin/surv_annot_process.R ${annot} ${summary} ${name_string}
		"""
	}

	process intersect_beds{
		label 'short_himem'
		label 'bedtools'
		input:
			val name_string from params.names
			file "ins.bed" from ch_ins
			file "inv.bed" from ch_inv
			file "del.bed" from ch_del
			file "dup.bed" from ch_dup
			file "tra.bed" from ch_tra
		output:
			file "${name_string}.ins.s.bed" into ch_ins_s
			file "${name_string}.ins.e.bed" into ch_ins_e
			file "${name_string}.del.s.bed" into ch_del_s
			file "${name_string}.del.s.bed" into ch_del_e
			file "${name_string}.inv.e.bed" into ch_inv_e
			file "${name_string}.tra.e.bed" into ch_tra_e
			file "${name_string}.dup.e.bed" into ch_dup_e
			file "${name_string}.ins.genes.bed" into ch_ins_genes
			file "${name_string}.del.genes.bed" into ch_del_genes
			file "${name_string}.inv.genes.bed" into ch_inv_genes
			file "${name_string}.dup.genes.bed" into ch_dup_genes
			file "${name_string}.tra.genes.bed" into ch_tra_genes
			file "${name_string}.ins.exons.bed" into ch_ins_exons
			file "${name_string}.del.exons.bed" into ch_del_exons
			file "${name_string}.inv.exons.bed" into ch_inv_exons
			file "${name_string}.dup.exons.bed" into ch_dup_exons
			file "${name_string}.tra.exons.bed" into ch_tra_exons
		script:
		"""
		/usr/bin/env bash ${projectDir}/bin/intersect_beds.sh ${name_string}
		"""
	}

	process summarize_intersections{
		publishDir params.outdir, mode:'copy'
		label 'short_himem'
		label 'tidyverse'
		input:
			val name_string from params.names
			file "summary.csv" from ch_summary
			file "ins.bed" from ch_ins
			file "dup.bed" from ch_dup
			file "tra.bed" from ch_tra
			file "inv.bed" from ch_inv
			file "del.bed" from ch_del
			file "ins.s.bed" from ch_ins_s
			file "ins.e.bed" from ch_ins_e
			file "del.s.bed" from ch_del_s
			file "del.e.bed" from ch_del_e
			file "dup.e.bed" from ch_dup_e
			file "inv.e.bed" from ch_inv_e
			file "tra.e.bed" from ch_tra_e
			file "ins.genes.bed" from ch_ins_genes
			file "del.genes.bed" from ch_del_genes
			file "dup.genes.bed" from ch_dup_genes
			file "tra.genes.bed" from ch_tra_genes
			file "inv.genes.bed" from ch_inv_genes
			file "ins.exons.bed" from ch_ins_exons
			file "del.exons.bed" from ch_del_exons
			file "dup.exons.bed" from ch_dup_exons
			file "tra.exons.bed" from ch_tra_exons
			file "inv.exons.bed" from ch_inv_exons
		output:
			file "${name_string}.survivor_results.csv"
		script:
		"""
		/usr/bin/env Rscript ${projectDir}/bin/post_filt.R
		mv survivor_results.csv ${name_string}.survivor_results.csv
		"""
	}
	
	process annotate_exons{
    publishDir params.outdir, mode:'copy'
    label 'midjob'
    label 'pysam'
    input:
        val name_string from params.names
        file "merged_vcf" from ch_merged_vcf
        file "ins.exons.bed" from ch_ins_exons
        file "del.exons.bed" from ch_del_exons
        file "dup.exons.bed" from ch_dup_exons
        file "tra.exons.bed" from ch_tra_exons
        file "inv.exons.bed" from ch_inv_exons
    output:
        file "${name_string}.merged.annotated.vcf"
    script:
    """
    /usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v merged_vcf \
        -i ins.exons.bed -d del.exons.bed \
        -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
        -o ${name_string}.merged.annotated.vcf
    """
	}
	*/

	// NEED TO CHANGE THIS TO HAVE CCS INCLUDED 
	
	
	//It is highly recommended to provide one tandem repeat annotation .bed file of your reference to pbsv discover via --tandem-repeats. This increases sensitivity and recall. Feel free to use the following for human SV calling: GRCh38 or hs37d5/hg19.
	
	// Call structural variants and assign genotypes
	// pbsv call ref.fa ref.sample1.svsig.gz ref.sample2.svsig.gz ref.var.vcf
	
	
}  // END OF if params.seqmode == 'pacbio'

// ======================================= Short Read Alignment Proccesses =====
if (params.seqmode == 'illumina') {
	if (params.fastq1){
	params.bwa = params.genome && params.fasta && params.fastq1 ? params.genomes[params.genome].bwa ?: null : null
	process BuildBWAindexes {
		tag "${fasta}"
		label 'bwa'
		publishDir params.outdir, mode: 'copy'
		input:
			file(fasta) from ch_fasta
		output:
			file("${fasta}.*") into bwa_built
		when: !(params.bwa) && params.fasta && params.fastq1
		script:
		"""
		bwa index ${fasta}
		"""
	}

	ch_bwa = params.bwa ? Channel.value(file(params.bwa)) : bwa_built
	process readgroup{
	  label 'python2'
	  input:
		file fq1 from ch_fastq1
		val sample_name from params.names		
	  output:
		file "${sample_name}.rg" into readgroup
	  script:
	  """
	  /usr/bin/env python ${projectDir}/bin/read_group_from_fastq.py -o ${sample_name}.rg $fq1
	  """
	}
 
	process map{
	  label 'bwa'
	  label 'cpus_8'
	  stageInMode 'copy'
	  input:
		file fq1 from ch_fastq1
		file fq2 from ch_fastq2
		file faidx from ch_bwa
		file fasta from ch_fasta
		file rgr from readgroup
		val sample_name from params.names
	  output:
		file "${sample_name}.sam" into ch_sam_map
	  when: !(params.bam) && params.fastq1
	  script:
	  """
	  bwa mem -K 100000000 -R \$(cat $rgr) -t ${task.cpus} -M ${fasta} $fq1 $fq2 > ${sample_name}.sam
	  """
	}

	process map2{
	  publishDir "${params.outdir}/alignments", mode:'copy'
	  label 'cpus_8'
	  label 'samtools'
	  input:
		file sam from ch_sam_map
		val sample_name from params.names
	  output:
		file "${sample_name}.bam" into ch_bam_map
	  script:
	  """
	  samtools sort --threads ${task.cpus} -m 30G $sam > ${sample_name}.bam
	  """
	}
	ch_bam_und = ch_bam_map 
	}else{
	ch_bam_und = params.bam? Channel.value(file(params.bam)) : null
	}
	// end of optional mapping steps. ch_bam_und (either from mapping or from input) is passed to post processing and calling

	process dedup{
	  publishDir "${params.outdir}/alignments", mode:'copy'
	  label 'gatk'
	  label 'cpus_8'
	  input:
		file bamf from ch_bam_und
		val sample_name from params.names
	  output:
	  	tuple sample_name, file("${bamf.baseName}.md.bam"), file("${bamf.baseName}.md.bam.bai") into ch_sample_name_bam_bai
		file "${bamf.baseName}.metrics" into bam_markd_m
	  script:
	  	markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
	  """
	  gatk --java-options ${markdup_java_options} \
			MarkDuplicates \
			--MAX_RECORDS_IN_RAM 50000 \
			--INPUT ${bamf} \
			--METRICS_FILE ${bamf.baseName}.metrics \
			--TMP_DIR . \
			--ASSUME_SORT_ORDER coordinate \
			--CREATE_INDEX true \
			--REMOVE_DUPLICATES true \
			--OUTPUT ${bamf.baseName}.md.bam
	
		mv ${bamf.baseName}.md.bai ${bamf.baseName}.md.bam.bai
	  """
	}
	
	process avoid_race_condition{
		stageInMode 'copy'
		label 'samtools'
		output:
  		file file into avoid_race_cond_1
		"""
    	echo 'Hello world!' > file
    	"""
	}

	process avoid_race_condition_2{
		stageInMode 'copy'
		label 'lumpy'
		input:
		file arc2 from avoid_race_cond_1
		output:
  		file file into avoid_race_cond_2
		"""
    	echo 'Hello world!' > file
    	"""
	}
	
	process avoid_race_condition_3{
		stageInMode 'copy'
		label 'breakdancer'
		input:
		file arc3 from avoid_race_cond_2
		output:
  		file file into avoid_race_cond_3
		"""
    	echo 'Hello world!' > file
    	"""
	}

	process avoid_race_condition_4{
		stageInMode 'copy'
		label 'manta'
		input:
		file arc4 from avoid_race_cond_3
		output:
  		file file into avoid_race_cond_4
		"""
    	echo 'Hello world!' > file
    	"""
	} 

	process avoid_race_condition_5{
		stageInMode 'copy'
		label 'cnmops'
		input:
		file arc5 from avoid_race_cond_4
		output:
  		file file into avoid_race_cond_5
		"""
    	echo 'Hello world!' > file
    	"""
	}

	process avoid_race_condition_6{
		stageInMode 'copy'
		label 'delly'
		input:
		file arc6 from avoid_race_cond_5
		output:
  		file file into avoid_race_cond_6
		"""
    	echo 'Hello world!' > file
    	"""
	}

	process bam_insertsize{
	  publishDir "${params.outdir}/alignments", mode:'copy'
	  label 'samtools'
	  input:
		tuple sample_name, bam_input, bam_index from ch_sample_name_bam_bai
	  output:
		file "insert_size.txt" into ins_size_ch
	  script:
	  """
	  samtools stats $bam_input |grep "^IS" |awk '{a = a + \$2*\$3; b = b + \$3}END{print int(a/b)}' > insert_size.txt
	  """
	}

	process fastaindex{
	  publishDir "${params.outdir}/alignments", mode:'copy'
	  label 'samtools'
	  input:
		file fasta from ch_fasta

	  output:
		file "${fasta}.fai" into fasta_fai_ch

	  script:
	  """
	  samtools faidx $fasta
	  """
	}

	ch_sample_name_bam_bai.into {
  		in_brkdncr
  		in_delly
  		in_lumpy
  		in_manta
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Lumpy SV ~~~~~

	process lumpy_mapping {
		tag "$sample_name"
		label 'lumpy'
		label 'cpus_8'
		publishDir "${params.outdir}/alignments/mapped_lumpy", pattern: "*_alignBWA_lumpy.bam", mode: 'copy'
		cpus params.threads

		input:
			tuple sample_name, bam_input, bam_index from in_lumpy

		output:
			tuple sample_name, path(bam_bwa_lumpy) into ( bam_bwa_lumpy_ch, bam_bwa_lumpy_splits_ch )
			tuple sample_name, path(dis_unsorted_bam) into dis_unsorted_bam_ch

		script:
			log.info "Mapping Lumpy (Map clipped reads, read group info, extract discordant alignments)"

			// Lumpy File Names:
			bam_name_sort      = sample_name + "_alignBWA_ReadNameSort"
			bam_name_sort_full = sample_name + "_alignBWA_ReadNameSort.bam"
			bam_bwa_lumpy      = sample_name + "_alignBWA_lumpy.bam"
			bam_bwa_lumpy_sort = sample_name + "_alignBWA_lumpySort.bam"
			dis_unsorted_bam   = sample_name + "_discordants.unsorted.bam"
			dis_sorted_bam     = sample_name + "_discordants.sorted.bam"
			split_unsorted_bam = sample_name + "_splitters.unsorted.bam"
			split_sorted_bam   = sample_name + "_splitters.sorted.bam"
			"""
			# Clipped_rc reads mapping to Genome
			samtools sort -n ${bam_input} -o ${bam_name_sort_full} -@ ${params.threads}
			# manual read group info
			samtools view -h ${bam_name_sort_full} \
			| samblaster --acceptDupMarks --excludeDups --addMateTags \
						--ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
			| samtools view -@ ${params.threads} -S -b - > ${bam_bwa_lumpy}
			# Extract the discordant pairedend alignments
			samtools view -@ ${params.threads} -b -F 1294 ${bam_bwa_lumpy} > ${dis_unsorted_bam}
		"""
	}

	//exit 0

	process lumpy_bwa_sort {
		tag "$sample_name"
		label 'picard'
		label 'cpus_8'

		input:
		tuple sample_name, path(bam_bwa_lumpy) from bam_bwa_lumpy_ch

		output:
		tuple sample_name, path("${bam_bwa_lumpy_sort}.ba[mi]") into bam_bwa_lumpy_sort_ch

		script:
		log.info "Mapping Lumpy (bam sort)"
		bam_bwa_lumpy_sort = sample_name + "_alignBWA_lumpySort" //+ .bam
		"""
		picard SortSam I=${bam_bwa_lumpy} O=${bam_bwa_lumpy_sort}.bam \
			SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		"""
	}
	process lumpy_discordant_sort {
		tag "$sample_name"
		label 'picard'
		label 'cpus_8'

		input:
		tuple sample_name, path(dis_unsorted_bam) from dis_unsorted_bam_ch

		output:
		tuple sample_name, path("${dis_sorted_bam}.ba[mi]") into dis_sorted_bam_ch

		script:
		log.info "Mapping Lumpy (discordant sort)"
		dis_sorted_bam = sample_name + "_discordants.sorted" //+ .bam
		"""
		picard SortSam I=${dis_unsorted_bam} O=${dis_sorted_bam}.bam \
			SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		"""
	}
	process lumpy_extract_splits {
		tag "$sample_name"
		label 'lumpy'
		label 'cpus_8'

		input:
		tuple sample_name, path(bam_bwa_lumpy) from bam_bwa_lumpy_splits_ch

		output:
		tuple sample_name, path("${split_unsorted_bam}.ba[mi]") into split_unsorted_bam_ch

		script:
		log.info "Mapping Lumpy (extract the splitread alignments)"
		extractSplitReads_BwaMem = "extractSplitReads_BwaMem"
		split_unsorted_bam = sample_name + "_splitters.unsorted" //+ .bam
		"""
		samtools view -h ${bam_bwa_lumpy} \
		| ${extractSplitReads_BwaMem} -i stdin \
		| samtools view -Sb - > ${split_unsorted_bam}.bam
		"""
	}
	process lumpy_split_bam_sort {
		tag "$sample_name"
		label 'picard'
		label 'cpus_8'
		publishDir "${params.outdir}/alignments/mapped_lumpy", pattern: "*_splitters.sorted.ba*", mode: 'copy'

		input:
		tuple sample_name, path(split_unsorted_bam) from split_unsorted_bam_ch

		output:
		tuple sample_name, path("${split_sorted_bam}.ba[mi]") into split_sorted_bam_ch

		script:
		log.info "Mapping Lumpy (sort split reads)"
		split_sorted_bam = sample_name + "_splitters.sorted" //+ .bam
		"""
		picard SortSam I=${split_unsorted_bam[0]} O=${split_sorted_bam}.bam \
			SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		"""
	}
	process lumpy_call_sv {
		tag "$sample_name"
		label 'lumpy'
		label 'cpus_8'
		publishDir "${params.outdir}/alignments/mapped_lumpy", pattern: "*_discordants.sorted.ba*", mode: 'move'

		//errorStrategy { task.exitStatus=141 ? 'ignore' : 'terminate' } // validExitStatus 141 for pairend_distro

		input:
			tuple sample_name, path(bam_bwa_lumpy_sort) from bam_bwa_lumpy_sort_ch
			tuple sample_name, path(split_sorted_bam) from split_sorted_bam_ch
			tuple sample_name, path(dis_sorted_bam) from dis_sorted_bam_ch
			val abs_outdir from abs_outdir

		output:
	//    path(dis_sorted_bam, includeInputs=true) 
	// NOTE: The above statement was causing issues in the run. Not entirely sure that we need the sorted bam as output. 
	//       It has been turned off for now. 
			file lumpy_sort_vcf into reheader_lumpy
			
			

		shell:
			log.info "Call SV by Lumpy, sort vcf"
			pairend_distro = "pairend_distro.py"
			histo          = sample_name + "_alignBWA_lumpySort.lib1.histo"
			lumpy_vcf      = sample_name + "_lumpyOut.vcf"
			lumpy_sort_vcf = sample_name + "_lumpySort.vcf"
			exclude_regions = "/ref_data/mm10.gaps.centro_telo.scafold.exclude.bed"
			'''
			RG_ID=$(samtools view -H !{bam_bwa_lumpy_sort[1]} | grep '^@RG' | sed "s/.*ID:\\([^\\t]*\\).*/\\1/g")
			#orig: metrics=$(samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>&1
			samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 > pre_metrics 2>/dev/null
			metrics=$(cat pre_metrics | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>/dev/null \
				&& [ $? = 141 ] && echo 'metrics to pairend_distro had exitcode: '$?;
			mean=$(echo "${metrics}" | cut -d " " -f 1)
			mean=$(echo "${mean}"    | cut -d ":" -f 2)
			std_dev=$(echo "${metrics}" | cut -d " " -f 2)
			std_dev=$(echo "${std_dev}" | cut -d ":" -f 2)
			rm pre_metrics;

			lumpy \
				-mw 4 \
				-x !{exclude_regions} \
				-pe id:"${RG_ID}",bam_file:!{dis_sorted_bam[1]},histo_file:!{histo},mean:"${mean}",stdev:"${std_dev}",read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
				-sr id:"${RG_ID}",bam_file:!{split_sorted_bam[1]},back_distance:10,weight:1,min_mapping_threshold:20 \
				> !{lumpy_vcf}

			vcfSort.sh !{lumpy_vcf} !{lumpy_sort_vcf}
			'''
	}

	process reheader_lumpy {
		tag "$sample_name"
		label 'bcftools'
		label 'tiny_job'
		publishDir "${params.outdir}/lumpySVOut", pattern: "*_lumpySort.vcf", mode: 'move'

		input:
			val abs_outdir from abs_outdir
			val sample_name from params.names
			file "lumpySort.vcf" from reheader_lumpy
		
		output:
			path("${sample_name}_lumpySort.vcf")
			path("vcf_path") into vcf_lumpy
		
		script:
			log.info "Reheading lumpy SV VCF"
			"""
			printf "${sample_name}_lumpy\n" > rehead_lumpy.txt
			bcftools reheader --samples rehead_lumpy.txt \
				-o ${sample_name}_lumpySort.vcf \
				lumpySort.vcf
			echo ${abs_outdir}/lumpySVOut/${sample_name}_lumpySort.vcf > vcf_path # for later merging
			"""
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BreakDancer SV ~~~~~
	process breakdancer_calling_sv {
    tag "$sample_name"
    label 'breakdancer'
	label 'cpus_8'

    input:
		tuple sample_name, bam_input, bam_index from in_brkdncr
		val abs_outdir from abs_outdir

    output:
		file breakdancer_sv_out into ch_breakdancer_sv

    script:
		breakdancer_config   = sample_name + "_config"
		breakdancer_sv_out   = sample_name + "_BreakDancer-SV"

		log.info "Calling BreakDancer SV"
		"""
		bam2cfg.pl ${bam_input} > ${breakdancer_config}
		breakdancer-max -r 5 -s 50 -h ${breakdancer_config} > ${breakdancer_sv_out}
		"""
	}
	process breakdancer_sv_to_vcf {
    tag "$sample_name"
    label 'python'
	label 'cpus_8'

    input:
		tuple sample_name, bam_input, bam_index from in_brkdncr
		val abs_outdir from abs_outdir
		file breakdancer_sv_out from ch_breakdancer_sv

    output:
		file breakdancer_sort_vcf into reheader_breakdancer

    script:
		breakdancer_2_vcf    = sample_name + "_BreakDancer2VCF.vcf"
		breakdancer_sort_vcf = sample_name + "_BreakDancerSortVCF.vcf"

		log.info "Converting BreakDancer SV calls to VCF"
		"""

		breakdancer2vcfHeader.py -i ${breakdancer_sv_out} -o ${breakdancer_2_vcf}
		vcfSort.sh ${breakdancer_2_vcf} ${breakdancer_sort_vcf}
		"""
	}

	process reheader_breakdancer {
		tag "$sample_name"
		label 'bcftools'
		label 'tiny_job'
		publishDir "${params.outdir}/BreakDancerSVOut", pattern: "*_BreakDancerSortVCF.vcf", mode: 'move'

		input:
			val abs_outdir from abs_outdir
			val sample_name from params.names
			file "BreakDancerSortVCF.vcf" from reheader_breakdancer
		
		output:
			path("${sample_name}_BreakDancerSortVCF.vcf")
			path("vcf_path") into vcf_breakdancer
		
		script:
			log.info "Reheading Breakdancer SV VCF"
			"""
			printf "${sample_name}_breakdancer\n" > rehead_breakdancer.txt
			bcftools reheader --samples rehead_breakdancer.txt \
				-o ${sample_name}_BreakDancerSortVCF.vcf \
				BreakDancerSortVCF.vcf
			echo ${abs_outdir}/BreakDancerSVOut/${sample_name}_BreakDancerSortVCF.vcf > vcf_path # for later merging
			"""
	}
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manta SV ~~~~~

	process manta_calling_sv {
		tag "$sample_name"
		label 'manta'
		label 'cpus_8'
		
		publishDir "${params.outdir}/temps", pattern: "mantaSVOut", enabled: params.keep_intermediate
		cpus params.threads

		input:
			tuple sample_name, bam_input, bam_index from in_manta
			file fasta from ch_fasta
			file fastfai from fasta_fai_ch
			val abs_outdir from abs_outdir

		output:
			file "candidateSV.vcf" into reheader_manta

		script:
			log.info "Calling Manta SV"
			"""
			/usr/local/bin/configManta.py \
				--runDir mantaSVOut \
				--bam ${bam_input} \
				--referenceFasta ${fasta}
			./mantaSVOut/runWorkflow.py -m local -j ${params.threads}
			mv mantaSVOut/results/variants/candidateSV.vcf.gz ./
			gunzip candidateSV.vcf.gz
			"""
	}

	process reheader_manta {
		tag "$sample_name"
		label 'bcftools'
		label 'tiny_job'
		publishDir "${params.outdir}/mantaSVout", pattern: "${sample_name}_candidateSV.vcf", mode: 'move'

		input:
			val abs_outdir from abs_outdir
			val sample_name from params.names
			file "candidateSV.vcf" from reheader_manta
		
		output:
			path("${sample_name}_candidateSV.vcf")
			path("vcf_path") into vcf_manta
		
		script:
			log.info "Reheading Manta SV VCF"
			"""
			printf "${sample_name}_manta\n" > rehead_manta.txt
			mv candidateSV.vcf ${sample_name}_candidateSV.vcf
			echo ${abs_outdir}/mantaSVout/${sample_name}_candidateSV.vcf > vcf_path # for later merging
			"""
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Delly SV ~~~~~

	process delly_calling_sv {
    tag "$sample_name"
    label 'delly'
	label 'cpus_8'

    input:
		tuple sample_name, bam_input, bam_index from in_delly
		file fasta from ch_fasta
		file fastfai from fasta_fai_ch
    
    output:
    	tuple sample_name, path(delly_bcf) into delly_bcf_out

    script:
    	delly_bcf = sample_name + "_Dellybcf.bcf"
    	log.info "Calling Delly SV"
		exclude_regions = "/ref_data/mouse.mm10.excl.tsv"
    """
    delly call \
          -q 40 \
          -x ${exclude_regions} \
          -s 500 \
          -o ${delly_bcf} \
          -g ${fasta} ${bam_input}
    """
	}
	process delly_bcf2vcf_sort {
		tag "$sample_name"
		label 'bcftools'
		label 'cpus_8'		

		input:
			tuple sample_name, path(delly_bcf) from delly_bcf_out
			val abs_outdir from abs_outdir

		output:
			file "${delly_sort_vcf}" into reheader_delly

		script:
			delly_vcf      = sample_name + "_DellyVCF.vcf"
			delly_sort_vcf = sample_name + "_dellySort.vcf"

			log.info "Delly bcf2vcf to Sorted VCF"
			"""
			bcftools view ${delly_bcf} > ${delly_vcf}
			vcfSort.sh ${delly_vcf} ${delly_sort_vcf}
			"""
	}

	process reheader_delly {
		tag "$sample_name"
		label 'bcftools'
		label 'tiny_job'
		publishDir "${params.outdir}/DellySVOut", pattern: "*_dellySort.vcf", mode: 'move'

		input:
			val abs_outdir from abs_outdir
			val sample_name from params.names
			file "dellySort.vcf" from reheader_delly
		
		output:
			path("${sample_name}_dellySort.vcf")
			path("vcf_path") into vcf_delly
		
		script:
			log.info "Reheading Delly SV VCF"
			"""
			printf "${sample_name}_delly\n" > rehead_delly.txt
			bcftools reheader --samples rehead_delly.txt \
				-o ${sample_name}_dellySort.vcf \
				dellySort.vcf
			echo ${abs_outdir}/DellySVOut/${sample_name}_dellySort.vcf > vcf_path # for later merging
			"""
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Collect all VCF file Paths  ~~~~~
	// NOTE: merged extension contains 'BDLM' => Breakdancer + Delly + Lumpy + Manta
	vcf_breakdancer.concat(
	vcf_delly,
	vcf_lumpy,
	vcf_manta
	)
	.collectFile(name: "sample.vcfs.txt", sort: false)
	.set { sample_vcfs_paths }
} // END OF if param.seqmode == 'illumina'
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Merge all SV calls VCF files  ~~~~~

process survivor{
	label 'long_himem'
	label 'survivor'
	input:
		path(vcf_paths) from sample_vcfs_paths
		val sample_name from params.names
		val surv_dist from params.surv_dist
		val surv_supp from params.surv_supp
		val surv_type from params.surv_type
		val surv_strand from params.surv_strand
		val surv_min from params.surv_min
	output:
		tuple sample_name, path("${sample_name}_mergedCall.BDLM.vcf") into ( vcf_merged, vcf_mrg_annot )
	script:
	"""
	SURVIVOR merge ${vcf_paths} ${surv_dist} ${surv_supp} ${surv_type} ${surv_strand} 0 ${surv_min} ${sample_name}_mergedCall.BDLM.vcf
	"""
}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Annotate if in Exon  ~~~~~

vcf_merged.into{vcf_merged_1; vcf_merged_2; vcf_merged_3}

process annotate_sv{
	label 'mid_job'
	publishDir params.outdir, mode:'copy'
	input:
		tuple sample_name, path(in_vcf) from vcf_merged_1
		val seqmode from params.seqmode
	output:
		file "${sample_name}.merged.overlap.annotated.txt" into ch_annot
	script:
	"""
	/usr/bin/env bash ${projectDir}/bin/surv_annot.sh ${sample_name} ${in_vcf} ${seqmode}
	"""
}


process summarize_sv{
	label 'mid_job'
	label 'pyvcf'
	input:
		tuple sample_name, path(in_vcf) from vcf_merged_2
	output:
		file "${sample_name}.survivor_summary.csv" into ch_summary
	script:
	"""
	/usr/bin/env python ${projectDir}/bin/sv_to_table.py -v ${in_vcf} -o ${sample_name}.survivor_summary.csv
	"""
}

ch_summary.into{ch_summary_1; ch_summary_2}

process prep_beds{
	label 'short_himem'
	label 'tidyverse'
	input:
		file annot from ch_annot
		file summary from ch_summary_1
		val name_string from params.names
	output:
		file "${name_string}.ins.bed" into ch_ins
		file "${name_string}.inv.bed" into ch_inv
		file "${name_string}.del.bed" into ch_del
		file "${name_string}.dup.bed" into ch_dup
		file "${name_string}.tra.bed" into ch_tra
	script:
	"""
	/usr/bin/env Rscript ${projectDir}/bin/surv_annot_process.R ${annot} ${summary} ${name_string}
	"""
}

ch_ins.into{ch_ins_1; ch_ins_2}
ch_inv.into{ch_inv_1; ch_inv_2}
ch_del.into{ch_del_1; ch_del_2}
ch_dup.into{ch_dup_1; ch_dup_2}
ch_tra.into{ch_tra_1; ch_tra_2}

process intersect_beds{
	label 'short_himem'
	label 'bedtools'
	input:
		val name_string from params.names
		file "ins.bed" from ch_ins_1
		file "inv.bed" from ch_inv_1
		file "del.bed" from ch_del_1
		file "dup.bed" from ch_dup_1
		file "tra.bed" from ch_tra_1
	output:
		file "${name_string}.ins.s.bed" into ch_ins_s
		file "${name_string}.ins.e.bed" into ch_ins_e
		file "${name_string}.del.s.bed" into ch_del_s
		file "${name_string}.del.s.bed" into ch_del_e
		file "${name_string}.inv.e.bed" into ch_inv_e
		file "${name_string}.tra.e.bed" into ch_tra_e
		file "${name_string}.dup.e.bed" into ch_dup_e
		file "${name_string}.ins.genes.bed" into ch_ins_genes
		file "${name_string}.del.genes.bed" into ch_del_genes
		file "${name_string}.inv.genes.bed" into ch_inv_genes
		file "${name_string}.dup.genes.bed" into ch_dup_genes
		file "${name_string}.tra.genes.bed" into ch_tra_genes
		file "${name_string}.ins.exons.bed" into ch_ins_exons
		file "${name_string}.del.exons.bed" into ch_del_exons
		file "${name_string}.inv.exons.bed" into ch_inv_exons
		file "${name_string}.dup.exons.bed" into ch_dup_exons
		file "${name_string}.tra.exons.bed" into ch_tra_exons
	script:
	"""
	/usr/bin/env bash ${projectDir}/bin/intersect_beds.sh ${name_string}
	"""
}

ch_ins_exons.into{ch_ins_exons_1; ch_ins_exons_2}
ch_del_exons.into{ch_del_exons_1; ch_del_exons_2}
ch_inv_exons.into{ch_inv_exons_1; ch_inv_exons_2}
ch_dup_exons.into{ch_dup_exons_1; ch_dup_exons_2}
ch_tra_exons.into{ch_tra_exons_1; ch_tra_exons_2}

process summarize_intersections{
	publishDir params.outdir, mode:'copy'
	label 'short_himem'
	label 'tidyverse'
	input:
		val name_string from params.names
		file "summary.csv" from ch_summary_2
		file "ins.bed" from ch_ins_2
		file "dup.bed" from ch_dup_2
		file "tra.bed" from ch_tra_2
		file "inv.bed" from ch_inv_2
		file "del.bed" from ch_del_2
		file "ins.s.bed" from ch_ins_s
		file "ins.e.bed" from ch_ins_e
		file "del.s.bed" from ch_del_s
		file "del.e.bed" from ch_del_e
		file "dup.e.bed" from ch_dup_e
		file "inv.e.bed" from ch_inv_e
		file "tra.e.bed" from ch_tra_e
		file "ins.genes.bed" from ch_ins_genes
		file "del.genes.bed" from ch_del_genes
		file "dup.genes.bed" from ch_dup_genes
		file "tra.genes.bed" from ch_tra_genes
		file "inv.genes.bed" from ch_inv_genes
		file "ins.exons.bed" from ch_ins_exons_1
		file "del.exons.bed" from ch_del_exons_1
		file "dup.exons.bed" from ch_dup_exons_1
		file "tra.exons.bed" from ch_tra_exons_1
		file "inv.exons.bed" from ch_inv_exons_1
	output:
		file "${name_string}.survivor_results.csv"
	script:
	"""
	/usr/bin/env Rscript ${projectDir}/bin/post_filt.R
	mv survivor_results.csv ${name_string}.survivor_results.csv
	"""
}
	
process annotate_exons{
publishDir params.outdir, mode:'copy'
label 'midjob'
label 'pysam'
input:
	tuple sample_name, path(in_vcf) from vcf_merged_3
	file "ins.exons.bed" from ch_ins_exons_2
	file "del.exons.bed" from ch_del_exons_2
	file "dup.exons.bed" from ch_dup_exons_2
	file "tra.exons.bed" from ch_tra_exons_2
	file "inv.exons.bed" from ch_inv_exons_2
output:
	file "${sample_name}.mergedCalls.InExon.BDLM.vcf"
script:
"""
/usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v ${in_vcf} \
	-i ins.exons.bed -d del.exons.bed \
	-u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
	-o ${sample_name}.mergedCalls.InExon.BDLM.vcf
"""
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Info ~ ~ ~ ~ ~ ~
/*
workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = workflow.complete
    wfEnd['Duration']     = workflow.duration
    wfEnd['Exit status']  = workflow.exitStatus
    wfEnd['Success']      = workflow.success
    if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = workflow.errorMessage
        wfEnd['.    Report']  = workflow.errorReport
    }
    Summary.show(wfEnd)
}
*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
