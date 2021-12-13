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
	  samtools sort --threads ${task.cpus} -m 2G ${sam} > ${name_string}.ngmlr.aligned.bam
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
		script:
			  """
			  sniffles -m ${bam} -v sniffles_calls.vcf
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
			"""
			pbsv call --ccs ${fasta} ${svsig} pbsv_calls.vcf
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
			"""
			pbsv call ${fasta} ${svsig} pbsv_calls.vcf
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

	process survivor{
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
	}

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


	// NEED TO CHANGE THIS TO HAVE CCS INCLUDED 
	
	
	//It is highly recommended to provide one tandem repeat annotation .bed file of your reference to pbsv discover via --tandem-repeats. This increases sensitivity and recall. Feel free to use the following for human SV calling: GRCh38 or hs37d5/hg19.
	
	// Call structural variants and assign genotypes
	// pbsv call ref.fa ref.sample1.svsig.gz ref.sample2.svsig.gz ref.var.vcf
	
	
}  // END OF if params.seqmode == 'pacbio'


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
	  output:
		file "${fq1.baseName}.rg" into readgroup
	  script:
	  """
	  /usr/bin/env python ${projectDir}/bin/read_group_from_fastq.py -o ${fq1.baseName}.rg $fq1
	  """
	}
 
	process map{
	  label 'bwa'
	  stageInMode 'copy'
	  input:
		file fq1 from ch_fastq1
		file fq2 from ch_fastq2
		file faidx from ch_bwa
		file fasta from ch_fasta
		file rgr from readgroup
	  output:
		file "${fq1.baseName}.sam" into ch_sam_map
	  when: !(params.bam) && params.fastq1
	  script:
	  """
	  bwa mem -K 100000000 -R \$(cat $rgr) -t ${task.cpus} -M ${fasta} $fq1 $fq2 > ${fq1.baseName}.sam
	  """
	}

	process map2{
	  publishDir params.outdir, mode:'copy'
	  label 'cpus_8'
	  label 'samtools'
	  input:
		file sam from ch_sam_map
	  output:
		file "${sam.baseName}.bam" into ch_bam_map
	  script:
	  """
	  samtools sort --threads ${task.cpus} -m 2G $sam > ${sam.baseName}.bam
	  """
	}
	ch_bam_und = ch_bam_map 
	}else{
	ch_bam_und = params.bam? Channel.value(file(params.bam)) : null
	}
	// end of optional mapping steps. ch_bam_und (either from mapping or from input) is passed to post processing and calling

	process dedup{
	  publishDir params.outdir, mode:'copy'
	  label 'gatk'
	  label 'cpus_8'
	  input:
		file bamf from ch_bam_und
	  output:
		file "${bamf.baseName}.md.bam" into ch_bam
		file "${bamf.baseName}.md.bam.bai" into ch_bam_bai
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
	  publishDir params.outdir, mode:'copy'
	  label 'samtools'
	  input:
		file bamf from ch_bam
	  output:
		file "insert_size.txt" into ins_size_ch
	  script:
	  """
	  samtools stats $bamf |grep "^IS" |awk '{a = a + \$2*\$3; b = b + \$3}END{print int(a/b)}' > insert_size.txt
	  """
	}

	process fastaindex{
	  publishDir params.outdir, mode:'copy'
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

	process lumpy{
	  publishDir params.outdir, mode:'copy'
	  label 'cpus_8'
	  label 'lumpy'
	  input:
		file bamf from ch_bam

	  output:
		file "lumpy_output.vcf" into lumpy_out
	  script:
	  """
	  lumpyexpress -B $bamf -P -m 2 -o lumpy_output.vcf -v
	  """
	}


	process breakdancer{
	  publishDir params.outdir, mode:'copy'
	  label 'cpus_8'
	  label 'breakdancer'
	  input:
		file bamf from ch_bam

	  output:
		file "${bamf.baseName}.cfg" into bd_cfg
		file "${bamf.baseName}.calls" into bd_calls
   
	  script:
	  """
	  /usr/bin/env perl bam2cfg.pl -q 30 -n 10000 $bamf > ${bamf.baseName}.cfg
	  breakdancer-max ${bamf.baseName}.cfg > ${bamf.baseName}.calls
	  """
	}


	process bd_to_vcf{
	  publishDir params.outdir, mode:'copy'
	  label 'cpu_1'
	  label 'python'

	  input:
		file calls from bd_calls

	  output:
		file "${calls.baseName}.vcf" into bd_vcf

	// This code was taken from the original SVE: https://github.com/TheJacksonLaboratory/SVE/blob/abb75e4a8269b75478e08ea18f05a2cb8b3b2521/stages/utils/breakdancer2vcf.py
	  script:
	  """
#!/usr/bin/env python
import time
from time import strftime

#input arguments: sys.argv[n]
#1-script.py, 2-refname 3-bd_calls.txt 4-outputname

# opens a valid file breakdancer output file
def read_breakdancer(path):
   file = open(path)  #open file connection
   h,raw,table = False,[],[] #init variables
   #scan file for header, read raw rows
   for line in file:
      if h:
         s = line.replace('\\n','')
         raw.append(s) 
      elif line.rfind('#Chr') != -1:
         s = line.replace('\\n','')
         h = True
   file.close()
   #split each row by the \\t and store in table variable
   for row in raw:
      table.append(row.split('\\t'))     
   #CTX type issue makes the number of columns 11 instead of 12....
#   for i in range(1,len(table)):
#      if (len(table[i-1]) != len(table[i])):
#         return "Break Dancer Format Error:\\nMissmatched Dimensions"
   return table

# writes a new .vcf file from a formatted table
def write_vcf(path,header,vcf_table):
    with open(path,'w') as f:
        s = header
        s = header + '\\n'.join(['\\t'.join(i) for i in vcf_table])
        f.write(s)        

#Build the VCF4.0 Header
def vcf_header(ref):
    form = '##fileformat=VCFv4.1\\n'
    date = '##fileDate=' + strftime('%Y%m%d',time.localtime()) + '\\n' 
    src  = '##source=BreakDancer_Max-1.4.5\\n'
    ref  = '##reference=' + ref + '\\n'
    inf1 = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\\n'
    inf2 = '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\\n'
    inf3 = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\\n'
    inf4 = '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n'
    alt1 = '##ALT=<ID=CNV,Description="Copy number variable region">\\n'
    alt2 = '##ALT=<ID=DEL,Description="Deletion">\\n'
    alt3 = '##ALT=<ID=INS,Description="Insertion">\\n'
    alt4 = '##ALT=<ID=TRA,Description="Translocation Event">\\n' #this should be updated to a BND event?
    alt5 = '##ALT=<ID=DUP,Description="Duplication Event">\\n'
    t_hd = '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'
    return (form+date+src+ref+inf1+inf2+inf3+inf4+alt1+alt2+alt3+alt4+alt5+t_hd)

#selects needed fields and converts what is needed to make a .vcf
#:::TO DO::: include additional CHR2 VCF tags in the INFO field
def build_vcf(table):
   vcf_table = []
   #main loop to work through the rows...
   for i in range(0,len(table)):
      CHR = table[i][0]
      POS = table[i][1]
      ID = 'breakdancer_' + str(i)
      REF = 'N'
      #convert breakdancer type to vcf type
      t = table[i][6]
      SVTYPE  = 'CNV' #default type
      if   t == 'DEL': SVTYPE = 'DEL'
      elif t == 'INS': SVTYPE = 'INS'
      elif t == 'INV': SVTYPE = 'INV'
      elif t == 'ITX': SVTYPE = 'DUP' #intrachromasomal translocation chrx->chrx
      elif t == 'CTX': SVTYPE = 'TRA' #interchromasomal translocation chrx->chry
      elif t == 'DUP': SVTYPE = 'DUP'
      elif t == 'CNV': SVTYPE = 'CNV'
      ALT = '<' + SVTYPE + '>'
      QUAL  = table[i][8]
      FILTER = 'PASS'
      END = table[i][4]
      SVLEN = table[i][7]
      INFO = 'END='+END+';SVTYPE='+SVTYPE+';SVLEN='+SVLEN+';IMPRECISE'
      vcf_table += [[CHR,POS,ID,REF,ALT,QUAL,FILTER,INFO]]
   max_seq = max([len(i[0]) for i in vcf_table])
   vcf_table = sorted(vcf_table,key=lambda x: (x[0].zfill(max_seq),int(x[1])))
   return vcf_table

#Test Code Here
#input arguments: sys.argv[n]
#0-script.py, 1-refname 2-bd_calls_dir
glob_path = '/Users/tbecker/Documents/CourseWork/16_2016_Spring/S4_clean/'
call = '$calls'
table = read_breakdancer(call)
write_vcf('${calls.baseName}.vcf',vcf_header('human_g1k_v37_decoy'),build_vcf(table))
	  """
	}
// END OF PYTHON SCRIPT

	process run_cnmops{
	  publishDir params.outdir, mode:'copy'
	  label 'cnmops'
	  label 'cpus_8'
	  stageInMode 'copy'

	  input:
		file bamf from ch_bam
		file bamfbai from ch_bam_bai
		file fastaf from ch_fasta

	  output:
		file "cnmops_output_calls.vcf" into cnmops_out
	  script:
	  """
	#!/usr/bin/env Rscript
	library(cn.mops)
	library(Rsamtools)
	mode <- 3
	in_bams <- strsplit("${bamf}", " ") [[1]]
	normal  <- 'quant'
	window <- 1000
	ref_min_len <- 10000 # Used to remove ERCC and other unlocalized
	cir_seg <- 'DNAcopy'
	out_vcf <- "cnmops_output_calls.vcf"
	outputdir = "./"
	outputprefix = "cnmops"

	f_names <- in_bams;
	cat('bam input files are: ',f_names,'\\n');

	s_names <- as.vector(matrix('',nrow=length(in_bams),ncol=1))
	for(i in 1:length(in_bams)){ 
		s_in_bam   <- strsplit(in_bams[i],'/')[[1]]
		s_in_bam   <- s_in_bam[length(s_in_bam)]
		s_names[i] <- strsplit(s_in_bam,'.bam')[[1]]
	}
	cat('sample labels are: ');
	cat(paste(s_names));

	#foreach chrom in the union of bam file set...
	ref <- c();
	for(i in 1:length(f_names)){
		header <- scanBamHeader(f_names[i]);
		ref <- union(ref,names(header[[1]]\$targets)[header[[1]]\$targets >= ref_min_len]);
	}

	res <- GRanges(); #empty place holder for error catching
	#calculate copy variation number regions   
	if(mode==0){ #mode = 0 -> cn.mops(regular multi sample WGS) 
		#try to make this more robust or alter input ref seq list = ref	
		tryCatch({
			cat('----------using cn.mops->multi sample WGS mode\\n');
			data <- getReadCountsFromBAM(BAMFiles=f_names, sampleNames=s_names, 
										 refSeqName = ref, WL=window);              
			cat(paste('----------finished reading ',dim(mcols(data))[1],' bins\\n',sep=''));
			#cleaning routine to suppress zero count issues.....cnmopsBUG
			for(contig in ref){ #remove from ref before runing anlaysis...
				if(sum(mcols(data[seqnames(data)%in%contig])[,1])==0){
					ref <- setdiff(ref,contig);
				}
			}
			cat('cleaned sequences are:')
			cat(paste(ref,sep=','))
			cat('\\n')
			data <- data[seqnames(data)%in%ref]; #adjust data using cleaned ref seqs
			seqlevels(data) <- ref;
			seqnames(data)  <- droplevels(seqnames(data));
			res  <- suppressWarnings(cn.mops(data,
									 I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
									 classes = paste('CN',seq(0,8,1),sep=''),
									 priorImpact = prior,
									 normType = normal, normQu = cutoff, norm = 1,
									 upperThreshold = upper, lowerThreshold = lower,
									 minWidth = min_seg, minReadCount = min_cnt)); #took out DNAcopy
			cat(paste('----------found ',dim(mcols(cnvr(res)))[1],' variations\n',sep=''));
		}, error = function(err){
			cat(paste(err));
			cat('\\n-------------------------------error----------------------------\\n')
			#save and debug the image here...
			image_m <- strsplit(out_vcf,'.vcf')[[1]]
			image_n <- sub('.bam','',image_m[length(image_m)])
			outimage <- paste(image_n,'.RData',sep='');
			save.image(file=outimage); #saves objects for remote debugging/ect..
		}, finally = {});
		cat(warnings());            	
	}

	if(mode==3){ #mode = 3 -> singlecn.mops(single sample WGS)
		#try to make this more robust or alter input ref seq list = ref	
			cat("BAM file: ", f_names, " ", s_names, " ", ref, " ", window)
			cat("\\n")
		tryCatch({
			cat('----------using singlecn.mops->single sample WGS mode\\n');
			data <- getReadCountsFromBAM(BAMFiles=f_names, sampleNames=s_names, 
										 refSeqName = ref, WL=window);              
			cat(paste('----------finished reading ',dim(mcols(data))[1],' bins\n',sep=''));
			#cleaning routine to suppress zero count issues.....cnmopsBUG
			for(contig in ref){ #remove from ref before runing anlaysis...
				if(sum(mcols(data[seqnames(data)%in%contig])[,1])==0){
					ref <- setdiff(ref,contig);
				}
			}
			cat('cleaned sequences are:')
			cat(paste(ref,sep=','))
			cat('\\n')
			data <- data[seqnames(data)%in%ref]; #adjust data using cleaned ref seqs
			seqlevels(data) <- ref;
			seqnames(data)  <- droplevels(seqnames(data));
			res  <- suppressWarnings(singlecn.mops(data,
									 I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
									 classes = paste('CN',seq(0,8,1),sep=''),
									 normType = normal, norm = 1))
			cat(paste('----------found ',dim(mcols(cnvr(res)))[1],' variations\\n',sep=''));
		}, error = function(err){
			cat(paste(err));
			cat('\\n-------------------------------error----------------------------\\n')
			#save and debug the image here...
			image_m <- strsplit(out_vcf,'.vcf')[[1]]
			image_n <- sub('.bam','',image_m[length(image_m)])
			outimage <- paste(image_n,'.RData',sep='');
			save.image(file=outimage); #saves objects for remote debugging/ect..
		}, finally = {});
		cat(warnings());         	
	}


	len <- 0;
	if(class(res)[1]=="CNVDetectionResult"){
		len <- dim(mcols(cnvr(res)))[1][1]; #number of CNV reported by cn.mops
	}
	if(len > 0){
		#estimate the integer copy number over the CNVRs and store in the CNVDetectionResult->res               
		res_icn <- calcIntegerCopyNumbers(res);                
		cnvrs <- cnvr(res_icn); #GRanges object
		chr <- as.vector(seqnames(cnvrs));
		cnr <- as.data.frame(ranges(cnvrs)); #get the start, end, width of CNV regions
		ics <- as.data.frame(mcols(cnvrs));
		icn <- as.data.frame(matrix(0,nrow=len,ncol=length(s_names)));
		colnames(ics) <- s_names;
		colnames(icn) <- colnames(ics);
		#convert 'CN0'->0, 'CN1'->1, etc...s_names<-labels has the correct info here
		for(i in 1:len){
			for(j in s_names){
				icn[i,j] <- as.numeric(strsplit(toString(ics[i,j]),'CN')[[1]][2]);
			}
		}
		ctype<- matrix('',nrow=len,ncol=1);
		colnames(icn) <- s_names; #sample names
		d_r <- dim(cnr)[1]; #number of cnv regions
		s_n <- dim(icn)[2]; #number of samples read
	
		#count up duplications, deletions or uncertain calls...IE genotyping...
		if(mode==3){ #each mode have to do del, dup or cnv
			dupn <- icn[,1] > 2;
			deln <- icn[,1] < 2;
			eqn  <- icn[,1] == 2;
		}
		if(mode==1){ #each mode have to do del, dup or cnv
			dupn <- icn[,1] < icn[,2];
			deln <- icn[,1] > icn[,2];
			eqn  <- icn[,1] == icn[,2];
		}
	
		for(i in 1:len){
			if(dupn[i]){ ctype[i] <- 'DUP'; }
			if(deln[i]){ ctype[i] <- 'DEL'; }
			if(eqn[i]) { ctype[i] <- 'CNV'; }
		}
		#convert 'CN0'->0, 'CN1'->1, etc...s_names<-labels has the correct info here
	
		#build the default vectors to be inserted into the final table
		CHROM <- matrix(chr,nrow=len, ncol=1);
		POS   <- cnr[,'start']; #copy over start values
		ID    <- matrix('',nrow=len, ncol=1);
		ID    <- as.matrix(paste('cnmops_',seq(1,len,1),ID,sep=''));
		REF   <- matrix('N',nrow=len, ncol=1); #N for a DEL, full seq ACGT...AAACT for ref
		ALT   <- matrix(paste('<',ctype,'>',sep=''),nrow=len, ncol=1);
		QUAL  <- matrix('.',nrow=len, ncol=1);
		FILTER<- matrix('PASS',nrow=len, ncol=1);
		INFO  <- matrix('IMPRECISE',nrow=len, ncol=1);
		FORMAT<- matrix('CN',nrow=len, ncol=1);
		#SAMPLES already have in result matrix
		SVTYPE<- gsub('<','',ALT);
		SVTYPE<- gsub('>','',SVTYPE);
		#make the INFO field with the meta tags and END position
		for(i in 1:len){
			INFO[i] <-paste('END=',cnr[i,'end'],';SVTYPE=',SVTYPE[i],';SVLEN=',
							cnr[i,'width'],';IMPRECISE',sep=''); #no seperation on values
		}

		#bind the colums into a new table
		table <- data.frame(cbind(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO),
							row.names=seq(1,len,1));
		d_t   <- dim(table)[2]; #default table size ~9
		#table[,(d_t+1):(d_t+s_n-1)] <- icn[2]; #insert integer copy numbers into the table
		#rename the colums for proper VCF compatibility
		c_n   <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'); 
		colnames(table) <- c_n; #attach the column names to the table
		#sort by chrom,pos the table if multi chrom...
		#table <- table[order(table[,'CHROM'],-table[,'POS']),];
		#redo the cmops_ids...
		#table[,'ID'] <- ID;
	} else { table <- data.frame(); } #empty data frame -> no calls...

	#VCF 4.0 Header construction-------------------------------------------------------
	# '.' => don't know the value like <NA> or unknown
	l01 <- '##fileformat=VCFv4.1\\n';
	l02 <- paste('##fileDate=',format(Sys.time(),'%Y%m%d'),'\\n',sep='');
	l03 <- '##source=cn.mops_CNV_calling\\n';
	l04 <- paste('##reference=',s_names[1],'\\n',sep='');
	l05 <- '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\\"Imprecise structural variation\">\\n';
	l06 <- '##INFO=<ID=END,Number=1,Type=Integer,Description=\\"End position of the variant described in this record\\">\\n';
	l07 <- '##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\\"Type of structural variant\\">\\n';
	l08 <- '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\\"Difference in length between REF and ALT alleles\\">\\n';
	l09 <- '##ALT=<ID=CNV,Description=\\"Copy number variable region\\">\\n';
	l10 <- '##ALT=<ID=DEL,Description=\\"Deletion\\">\\n';
	l11 <- '##ALT=<ID=DUP,Description=\\"Duplication\\">'; #put back \\n if adding format
	#l12 <- '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">';
	header <- paste(l01,l02,l03,l04,l05,l06,l07,l08,l09,l10,l11,sep='');
	#VCF 4.0 Header construction-------------------------------------------------------

	#write the header and then append the finished table
	#outvcf <- paste(outputdir,outputprefix,'.vcf',sep='');
	write(header,file=out_vcf,append=F);
	write.table(table,file=out_vcf,quote=F,append=T,col.names=T,row.names=F,sep='\\t');

	#save the plots...
	#seg_plot <- paste(outputdir,outputprefix,'_segplot.pdf',sep='');
	#pdf.options(width=8, height=10.5,onefile=T,paper='letter')
	#pdf(file=seg_plot)
	#segplot(res, seqname = chrom, sampleIdx = (1:n_sample))
	#dev.off()

	#read_plot <- paste(outputdir,outputprefix,'_readplot.pdf',sep='');	
	#plot(res, which = (1:n_sample),toFile=T,filename=read_plot)

	cat('\\ncn.mops cnv calling complete........')
	#END of cn.mops() cnv calling and visualization ===================================

	"""
	}
	// END OF RSCRIPT

	process manta_calling_sv {
		label 'manta'
		label 'cpus_8'
		
		publishDir params.outdir, mode: 'copy'

		input:
		file bamf from ch_bam
		file bamfb from ch_bam_bai
		file fasta from ch_fasta
		file fastafai from fasta_fai_ch		
		
		output:
		path("*candidateSV.vcf*")

		script:
		log.info "Calling Manta SV"
		"""
		/usr/bin/env python configManta.py \
			--runDir mantaSVOut \
			--bam ${bamf} \
			--referenceFasta ${fasta}
		/usr/bin/env python ./mantaSVOut/runWorkflow.py -m local -j ${task.cpus}
		mv mantaSVOut/results/variants/candidateSV.vcf.gz ./manta_candidateSV.vcf.gz
		"""
	}

	process hydra{
	  label 'hydra'
	  label 'cpus_8'
	  stageInMode 'copy'
	  input:
		file bamf from ch_bam
		file bamif from ch_bam_bai
	  output:
		path "routed/*" into hydra_out_all
		file "samples.conf" into conf_out, conf_out2
		file "$bamf" into bamc, bamc2
		file "$bamif" into bamci, bamci2
		file "${bamf}.bedpe" into bampe, bampe2
	  script:
	  """
	  echo "sample1\t$bamf" > samples.txt
	  /usr/bin/env python make_hydra_config.py -i samples.txt -s 100000 -n 16 > samples.conf
	  /usr/bin/env python extract_discordants.py -c samples.conf -d sample1 
	  mkdir routed
	  cd routed
	  ln -s ../${bamf} .
	  ln -s ../${bamif} .
	  ln -s ../${bamf}.bedpe .
	  hydra-router -config ../samples.conf -routedList ../samples.routed
	  rm ${bamf}
	  rm ${bamif}
	  rm ${bamf}.bedpe
	  #assemble-routed-files.sh samples.conf samples.routed 1 60
	  """
	}

	process hydra_asm{
	  label 'hydra'
	  label 'cpu_1'
	  input:
		each file(rt) from hydra_out_all
		file conf from conf_out
		file bamf from bamc
		file bamif from bamci
		file bamp from bampe
	  output:
		file "${rt}.posSorted.posClusters.assembled" into asm_out   
		file "${rt}.posSorted.posClusters.punted" into pun_out   

	  script:
	  """
	  echo "$rt"
	  hydra-assembler -config $conf -routed "$rt" -maxMappings 60
	  """
	}

	process hydra_collect{
	  publishDir params.outdir, mode:'copy'
	  label 'hydra'
	  label 'cpus_8'
	  input:
		file asm from asm_out.collect()
		file punt from pun_out.collect()
		file conf from conf_out2
		file bamf from bamc2
		file bamif from bamci2
		file bamp from bampe2

	  output:
		file "all.assembled" into asm_all
		file "hydra.calls.final" into asm_final
		file "hydra*" into asm_sv
	  script:
	  """
	  /usr/bin/env bash combine-assembled-files.sh . all.assembled
	  /usr/bin/env python forceOneClusterPerPairMem.py -i all.assembled -o hydra.calls
	  /usr/bin/env python frequency.py -c $conf -f hydra.calls.final -d hydra.calls.detail > hydra.calls.freq
	  grep -v "#" all-sv.calls.freq | /usr/bin/env python hydraToBreakpoint.py -i stdin > hydra.calls.bkpts
	  """
	}

	process fatobit{
	  label 'ucsc'
	  input:
		file fasta from ch_fasta  
	  output:
		file "${fasta}.2bit" into ch_2bit
	  script:
	  """
	  faToTwoBit $fasta ${fasta}.2bit
	  """
	}

	process hydra_vcf{
	  publishDir params.outdir, mode:'copy'
	  label 'bxpython'  
	  input:
		file fasta from ch_2bit
		file asm from asm_final
	  output:
		file "${asm.baseName}.vcf" into hydra_vcf
	  script:
	  """
	  /usr/bin/env python hydra_to_vcf.py $asm $fasta
	  """
	}

	// Run delly
	params.delly_exc = params.delly_genome ? "${params.delly_base}${params.delly_genome}.excl.tsv" : null
	dex = params.delly_exc ? Channel.value(file(params.delly_exc)) : "null"
	//type_list = Channel.from('DEL', 'DUP', 'INV', 'BND', 'INS')
	process delly{
	  label 'cpus_8'
	  label 'delly'
	  input:
	 //   val type from type_list
		file fasta from ch_fasta
		file bamf from ch_bam
		file bami from ch_bam_bai
		file excl from dex
	  output:
		file "delly.bcf" into delly_bcf
	  script:
	  exc = params.delly_exc ? "-x ${excl}" : ""
	  """
	  delly call -g $fasta -o delly.bcf $bamf 
	  """
	}

	process bcf2vcf{
	  publishDir params.outdir, mode:'copy'
	  label 'bcftools'
	  label 'cpus_8'
	  input:
		file bcf from delly_bcf
	  output:
		file "${bcf.baseName}.vcf" into delly_vcf
	  script:
	  """
	  bcftools view delly.bcf > delly.vcf
	  """
	}

	process cnvpytor{
	  publishDir params.outdir, mode:'copy'
	  label 'cnvpytor'
	  label 'cpus_8'
	  input:
		file bamf from ch_bam
		file bami from ch_bam_bai
		file fasta from ch_fasta
	  output:
		file "cnvpytor.calls" into cnvpytor_calls
	  script:
	  """
	  #cnvnator -root tree.root -tree $bamf
	  #cnvnator -his ${params.cnvpytor_bin} -root tree.root -fasta $fasta
	  #cnvnator -root tree.root -stat ${params.cnvpytor_bin}
	  #cnvnator -root tree.root -partition ${params.cnvpytor_bin}
	  #cnvnator -root tree.root -call ${params.cnvpytor_bin} > cnvnator.calls
	  
	  cnvpytor -root file.pytor -rd $bamf
	  cnvpytor -root file.pytor -his ${params.cnvpytor_bin}
	  cnvpytor -root file.pytor -partition ${params.cnvpytor_bin}
	  cnvpytor -root file.pytor -call ${params.cnvpytor_bin} > cnvpytor.calls
	  
	  
	  """
	}

	process cnv2vcf{
	  publishDir params.outdir, mode:'copy'
	  label 'perl'
	  label 'high_mem'
	  input:
		file calls from cnvpytor_calls
		file fasta from ch_fasta
	  output:
		file "cnvnator.vcf"
	  script:
	  """
	  /usr/bin/env perl cnvnator2VCF.pl -reference $fasta $calls ${params.genome} > cnvnator.vcf
	  
	  """
	}
} // END OF if param.seqmode == 'illumina'
