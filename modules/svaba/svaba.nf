process SVABA {
  tag "$sampleID"

  cpus = 8
  memory { normal_bam.size() < 60.GB ? 15.GB : 48.GB }
  time { normal_bam.size() < 60.GB ? '10:00:00' : '24:00:00' }
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/svaba:v0.2.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$sampleID" + '/callers'  : 'svaba' }", pattern: "*.vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

  output:
  tuple val(sampleID), path("*svaba.germline.indel.vcf.gz"), val(meta), val(normal_name), val(tumor_name), val('svaba'), emit: svaba_germline_indel_vcf
  tuple val(sampleID), path("*svaba.germline.sv.vcf.gz"), val(meta), val(normal_name), val(tumor_name), val('svaba'), emit: svaba_germline_sv_vcf
  tuple val(sampleID), path("*svaba.somatic.indel.vcf.gz"), val(meta), val(normal_name), val(tumor_name), val('svaba'), emit: svaba_somatic_indel_vcf
  tuple val(sampleID), path("*svaba.somatic.sv.vcf.gz"), val(meta), val(normal_name), val(tumor_name), val('svaba'), emit: svaba_somatic_sv_vcf
  tuple val(sampleID), path("*svaba.bps.txt.gz"), val(meta), val(normal_name), val(tumor_name), val('svaba'), emit: svaba_unfiltered_variants
  tuple val(sampleID), path("*svaba.contigs.bam"), emit: svaba_contigs_bam
  tuple val(sampleID), path("*svaba.discordant.txt.gz"), emit: svaba_discordants
  tuple val(sampleID), path("*svaba.log"), emit: svaba_log
  tuple val(sampleID), path("*svaba.alignments.txt.gz"), emit: svaba_alignments

  script:
  """
  svaba run \
    -t ${tumor_bam} \
    -n ${normal_bam} \
    -p ${task.cpus} \
    -a ${sampleID}_svaba \
    -G ${params.combined_reference_set} \
    --region ${params.callRegions} \
    -D ${params.dbSNP} \
    -z on
  """
}
// NOTE: VCF Output header has the BAM file names as 'sampleID' e.g.,: 
// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	test-test_realigned_BQSR.bam	test-test2_realigned_BQSR.bam

// Usage: svaba run -t <BAM/SAM/CRAM> -G <reference> -a myid [OPTIONS]

//   Description: SV and indel detection using rolling SGA assembly and BWA-MEM realignment

//   General options
//   -v, --verbose                        Select verbosity level (0-4). Default: 0
//   -h, --help                           Display this help and exit
//   -p, --threads                        Use NUM threads to run svaba. Default: 1
//   -a, --id-string                      String specifying the analysis ID to be used as part of ID common.
//   Main input
//   -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM.
//   -t, --case-bam                       Case BAM/CRAM/SAM file (eg tumor). Can input multiple.
//   -n, --control-bam                    (optional) Control BAM/CRAM/SAM file (eg normal). Can input multiple.
//   -k, --region                         Run on targeted intervals. Accepts BED file or Samtools-style string
//       --germline                       Sets recommended settings for case-only analysis (eg germline). (-I, -L5, assembles NM >= 3 reads)
//   Variant filtering and classification
//       --lod                            LOD cutoff to classify indel as non-REF (tests AF=0 vs AF=MaxLikelihood(AF)) [8]
//       --lod-dbsnp                      LOD cutoff to classify indel as non-REF (tests AF=0 vs AF=MaxLikelihood(AF)) at DBSnp indel site [5]
//       --lod-somatic                    LOD cutoff to classify indel as somatic (tests AF=0 in normal vs AF=ML(0.5)) [2.5]
//       --lod-somatic-dbsnp              LOD cutoff to classify indel as somatic (tests AF=0 in normal vs AF=ML(0.5)) at DBSnp indel site [4]
//       --scale-errors                   Scale the priors that a site is artifact at given repeat count. 0 means assume low (const) error rate [1]
//   Additional options
//   -L, --mate-lookup-min                Minimum number of somatic reads required to attempt mate-region lookup [3]
//   -s, --disc-sd-cutoff                 Number of standard deviations of calculated insert-size distribution to consider discordant. [3.92]
//   -c, --chunk-size                     Size of a local assembly window (in bp). Set 0 for whole-BAM in one assembly. [25000]
//   -x, --max-reads                      Max total read count to read in from assembly region. Set 0 to turn off. [50000]
//   -C, --max-coverage                   Max read coverage to send to assembler (per BAM). Subsample reads if exceeded. [500]
//       --no-interchrom-lookup           Skip mate lookup for inter-chr candidate events. Reduces power for translocations but less I/O.
//       --discordant-only                Only run the discordant read clustering module, skip assembly.
//       --num-assembly-rounds            Run assembler multiple times. > 1 will bootstrap the assembly. [2]
//       --num-to-sample                  When learning about inputs, number of reads to sample. [2,000,000]
//       --hp                             Highly parallel. Don't write output until completely done. More memory, but avoids all thread-locks.
//   Output options
//   -z, --g-zip                          Gzip and tabix the output VCF files. [off]
//   -A, --all-contigs                    Output all contigs that were assembled, regardless of mapping or length. [off]
//       --read-tracking                  Track supporting reads by qname. Increases file sizes. [off]
//       --write-extracted-reads          For the case BAM, write reads sent to assembly to a BAM file. [off]
//   Optional external database
//   -D, --dbsnp-vcf                      DBsnp database (VCF) to compare indels against
//   -B, --blacklist                      BED-file with blacklisted regions to not extract any reads from.
//   -Y, --microbial-genome               Path to indexed reference genome of microbial sequences to be used by BWA-MEM to filter reads.
//   -V, --germline-sv-database           BED file containing sites of known germline SVs. Used as additional filter for somatic SV detection
//   -R, --simple-seq-database            BED file containing sites of simple DNA that can confuse the contig re-alignment.
//   Assembly and EC params
//   -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength
//   -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA. 0 is fast, but requires error correcting. [0]
//   -K, --ec-correct-type                (f) Fermi-kit BFC correction, (s) Kmer-correction from SGA, (0) no correction (then suggest non-zero -e) [f]
//   -E, --ec-subsample                   Learn from fraction of non-weird reads during error-correction. Lower number = faster compute [0.5]
//       --write-asqg                     Output an ASQG graph file for each assembly window.
//   BWA-MEM alignment params
//       --bwa-match-score                Set the BWA-MEM match score. BWA-MEM -A [2]
//       --gap-open-penalty               Set the BWA-MEM gap open penalty for contig to genome alignments. BWA-MEM -O [32]
//       --gap-extension-penalty          Set the BWA-MEM gap extension penalty for contig to genome alignments. BWA-MEM -E [1]
//       --mismatch-penalty               Set the BWA-MEM mismatch penalty for contig to genome alignments. BWA-MEM -b [18]
//       --bandwidth                      Set the BWA-MEM SW alignment bandwidth for contig to genome alignments. BWA-MEM -w [1000]
//       --z-dropoff                      Set the BWA-MEM SW alignment Z-dropoff for contig to genome alignments. BWA-MEM -d [100]
//       --reseed-trigger                 Set the BWA-MEM reseed trigger for reseeding mems for contig to genome alignments. BWA-MEM -r [1.5]
//       --penalty-clip-3                 Set the BWA-MEM penalty for 3' clipping. [5]
//       --penalty-clip-5                 Set the BWA-MEM penalty for 5' clipping. [5]
