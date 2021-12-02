# Mouse Mutant Resource Structural Variant Detection

This repository holds a [Nextflow](https://www.nextflow.io/) implementation of structural variant callers for both short-read Illumina based data, and PacBio CCS and CLR data. It was developed to provide results that are compatible with the [Mouse Mutant Resource Database](https://mmrdb.jax.org).

For short-read data, either raw- or mapped reads can be used as input. 

If provided aligned short-reads, the alignment file is passed to the individual tools for SV calling. If supplied raw short-read data, raw-reads are mapped to a provided reference with [BWA mem](http://bio-bwa.sourceforge.net/bwa.shtml), duplicates are marked with GATK MarkDuplicates and the subsequent marked aligned BAM file is passed to the individual tools for SV calling. The short-reads SV callers included are: [cn.MOPS](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html), [Lumpy](https://github.com/arq5x/lumpy-sv), [CNVpytor](https://github.com/abyzovlab/CNVpytor) (the python version of [CNVnator](https://github.com/abyzovlab/CNVnator)), [Breakdancer](https://github.com/genome/breakdancer), [Hydra](https://github.com/arq5x/Hydra), [Delly](https://github.com/dellytools/delly), and [Manta](https://github.com/Illumina/manta).

For PacBio data, raw-reads are used as input. 

The PacBio based long-read data are mapped to a provided reference with [PBMM2](https://github.com/PacificBiosciences/pbmm2) and [NGMLR](https://github.com/philres/ngmlr). Mapped reads from PBMM2 are passed to [PBSV](https://github.com/PacificBiosciences/pbsv) for SV calling, and mapped reads from NGMLR are passed to [CuteSV](https://github.com/tjiangHIT/cuteSV), [Sniffles](https://github.com/fritzsedlazeck/Sniffles), and [SVIM](https://github.com/eldariont/svim).

For PacBio data, called SVs from `PBSV` and `Sniffles` are merged with `SURVIVOR`, and compared with mm10 SV annotations from Ensembl and Sanger to generate a table that shows whether each identified SV is present in previously published results and/or annotated genes and exons.

All tools for both short- and long-read data are obtained from publicly available software container repositories (e.g., [quay.io/biocontainers](quay.io/biocontainers)). We have included the ability to use a local FASTA based genome reference file, or to pull available genome references from iGenomes ([https://support.illumina.com/sequencing/sequencing_software/igenome.html](https://support.illumina.com/sequencing/sequencing_software/igenome.html)). 

## Running mmrSVD

### [Install Nextflow](https://www.nextflow.io/index.html#GetStarted)

Nextflow runs as a standalone application, which can be downloaded with the following command: 

```
curl -s https://get.nextflow.io | bash
```

Once you have downloaded Nextflow, ensure it that the tool is located in directory in your PATH, or when calling `Nextflow` ensure you use the full path to the executable. For example `/home/$USER/nextflow` 

### Install mmrSVD

```
git clone https://github.com/TheJacksonLaboratory/mmrSVD.git
```

### Run mmrSVD

Prior to starting the workflow, you must create a parameter file, or modified the example provided. The parameter file requires you specify `seqmode = "illumina"` or `seqmode = "pacbio"` if seqmode is pacbio, you must specify `clr` or `css` as the data type.  
	

For `pacbio` data, `fastq1` is used as the param for inputting the raw-reads, and `fastq2` is not used. 

`fasta` is used to specify reference data. You can provide an absolute path to a reference fasta, or use `--genome <ref_id>` in the workflow submission statement (see below).  

#### Parameter File Example

Below is a minimal example of the required parameters to run this pipeline using PacBio CCS data. A larger example that includes many optional parameters [can be found here](https://github.com/TheJacksonLaboratory/mmrsvd/blob/main/params.config).

```
params {
  // REQUIRED PARAMETERS -------------------------------------------------------
  fastq1 = "path_to_fastq1_input"
  outdir = "path_to_output_files"
  seqmode = "pacbio"
  pbmode = "ccs"
  names = "name_string"
}
``` 

#### Example Slurm Launch Wrapper: 

The workflow itself is designed for use in `SLURM` based HPC environments with `Singularity`, or in local environments with `Singularity` availability. See [https://singularity.lbl.gov/](https://singularity.lbl.gov/) for Singularity install information. 

The HPC Sumner has the required modules for running the workflow. The workflow can be launched using a wrapper script e.g., `run_mmrSVD_slurm.sh` 

```
#!/bin/bash
#SBATCH --job-name=nf_mmrSVD
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 7-00:00:00
#SBATCH --mem=2000

export NXF_SINGULARITY_CACHEDIR=/path/to/sing_cache
	
module load singularity
	
nextflow -c /path/to/params.config run /path/to/mmrSVD/main.nf -profile slurm,singularity --genome mm10
```
 
**NOTE 1: We are specifying a cache location for singularity images. This is optional**
`export NXF_SINGULARITY_CACHEDIR=/path/to/sing_cache`

**NOTE 2: We are specifying a genome version to use from igenomes with `--genome mm10`. For a list of available genomes is available in the [igenomes configuration file](https://github.com/TheJacksonLaboratory/mmrsvd/blob/main/conf/igenomes.config)**

**NOTE 3: To run the above example, `nextflow` must be in your path. Otherwise, specify the absolute path of the executable**

**NOTE 4: You must specify the paths to `params.config` and the location of `main.nf`**

## Long read output files

### Alignments
- ${name_string}.fastq.ngmlr.aligned.bam - BAM output from `ngmlr`
- ${name_string}.fastq.ngmlr.aligned_svim/ - Directory of alignment files from `SVIM`
- ${name_string}.fastq.pbmm2.aligned.bam - BAM alignment from `pbmm2` (part of `PBSV` workflow)

### Structural Variant Calls

- cutesv_calls.vcf - SV calls from `cutesv`
- pbsv_calls.vcf - SV calls from `PBSV`
- sniffles_calls.vcf - SV calls from `Sniffles`
- svim_variants.vcf - SV calls from `SVIM`
- ${name_string}.merged.annotated.vcf - Structural Variants from `PBSV` and `sniffles` merged with `SURVIVOR`. An additional INFO field `InExon=TRUE/FALSE` has been added to reflect whether the SV is contained within exonic regions of the mm10 genome.

### Summary table
- ${name_string}.survivor_results.csv - Summary table of the `SURVIVOR` merged VCF that lists the chromosome, position, SV name used by `SURVIVOR`, the SV type (insertion, deletion, duplication, inversion, translocation), whether the SV is present in Sanger annotations for mm10 (Y/N), whether the SV is present in Ensembl annotations for mm10 (Y/N), and the name of the gene and/or exons that the SV overlaps.