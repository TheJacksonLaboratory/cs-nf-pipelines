## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/deconstructsigs@sha256:009ddb6ed3ec2a0290a88b1e7027dd3caac1a2f5f3df3e8f68f410481d9323a3
## It is reproduced below as it exists there without modification

## Run deconstructSigs on a single sample's within a VCF
library(deconstructSigs)
library(optparse)
library(VariantAnnotation)



## Parse arguments
option_list = list(
  make_option(c("-f","--file"),       type="character", default=NULL,     help="[REQUIRED] Input VCF", metavar='vcf'),
  make_option(c("-r","--ref"),        type="character", default="GRCh38", help="[REQUIRED] Reference genome: [GRCh37,GRCh38]", metavar="character"),
  make_option(c("-c","--cosmic"),     type="character", default=NULL,     help="[REQUIRED] COSMIC reference file (.rda format)", metavar="rda"),
  make_option(c("-x","--highconf"),   type="logical",   default=FALSE,    help="When TRUE, consider high confidence mutations. Requires 'HighConfidence' flag in VCF INFO field.", metavar='logical'),
  make_option(c("-o","--output"),     type="character", default=NULL,     help="[REQUIRED] output prefix", metavar="character"),
  make_option(c("-s","--samplename"), type="character", default=NULL,     help="[REQUIRED] sample name. Shoud match the sample name in the input VCF", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Check arguments (for cases where we're using in a non-WDL context)
if (is.null(opt$file))                  stop("--file option is required")
if (is.null(opt$ref))                   stop("--ref option is required")
if (!opt$ref %in% c('GRCh37','GRCh38')) stop("Allowable values for --ref are: GRCh37,GRCh38 ")
if (is.null(opt$output))                stop("--output option is required")
if (is.null(opt$samplename))            stop("--samplename option is required")



###############
## Read data ##
###############

## Read cosmic reference data
load(opt$cosmic)


## Load appropriate reference genome 
if(opt$ref == "GRCh38"){
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  sig_ref = signatures.genome.cosmic.v3.2.march2021.grch38
  bsg_ref = BSgenome.Hsapiens.UCSC.hg38
  
} else if (opt$ref =="GRCh37"){
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  sig_ref = signatures.genome.cosmic.v3.2.march2021.grch37
  bsg_ref = BSgenome.Hsapiens.UCSC.hg19
  
} 


## Read input VCF, optionally filter for high confidence mutations 
vcf = readVcf(opt$file)

if(opt$highconf){
  vcf = vcf[info(vcf)$HighConfidence]
}



###########################################
## Prepare deconstructSigs input and run ##
###########################################

sbs96.input = data.frame()
tmb = c()

input = data.frame(seqnames(vcf),start(vcf),ref(vcf),alt(vcf))
input[,"SAMPLE"] = opt$samplename

sigs.input = mut.to.sigs.input(input, 
                                sample.id="SAMPLE", 
                                chr="seqnames.vcf.", 
                                pos="start.vcf.", 
                                ref="ref.vcf.", 
                                alt="value",
                                bsg=bsg_ref,
                                sig.type='SBS')
sbs96.input = rbind(sbs96.input, sigs.input)
tmb = append(tmb, nrow(input))
  

## Run deconstructSigs
sbs96.output = data.frame()
sbs96.tumor = data.frame()
sbs96.product = data.frame()
sbs96.diff = data.frame()

for(samplename in rownames(sbs96.input)){
  
  sbssigs = whichSignatures(tumor.ref=sbs96.input, 
                            signatures.ref=sig_ref, 
                            sample.id=samplename, 
                            contexts.needed=TRUE, 
                            signature.cutoff=0.00, 
                            tri.counts.method='default')
  
  sbs96.output = rbind(sbs96.output, sbssigs$weights)
  sbs96.tumor = rbind(sbs96.tumor, sbssigs$tumor)
  sbs96.product = rbind(sbs96.product, sbssigs$product)
  sbs96.diff = rbind(sbs96.diff, sbssigs$diff)
  
}

## Compute absolute mutation count assigned to each signature
sbs96.tmb.output = sweep(sbs96.output, 1, tmb, FUN="*")



############
## Output ##
############

## Prepare output
sbs96.output = cbind(Sample = rownames(sbs96.output), sbs96.output)
sbs96.tmb.output = cbind(Sample = rownames(sbs96.tmb.output), sbs96.tmb.output)
sbs96.tumor = cbind(Sample = rownames(sbs96.tumor), sbs96.tumor)
sbs96.product = cbind(Sample = rownames(sbs96.product), sbs96.product)
sbs96.diff = cbind(Sample = rownames(sbs96.diff), sbs96.diff)

## Write output 
write.table(sbs96.output,     file=paste0(opt$output,".txt"),               sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.tmb.output, file=paste0(opt$output,".counts.txt"),        sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.tumor,      file=paste0(opt$output,".input.txt"),         sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.product,    file=paste0(opt$output,".reconstructed.txt"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.diff,       file=paste0(opt$output,".diff.txt"),          sep="\t", quote=FALSE, row.names=FALSE)
