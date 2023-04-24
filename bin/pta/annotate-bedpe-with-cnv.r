## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Annotate a merged & annotated BEDPE with closest CNV changepoint
libs = c('optparse', 'StructuralVariantAnnotation', 'VariantAnnotation', 'rtracklayer', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)



## Handle non-standard bedpe columns better
readBEDPE = function(f) {
  
  ## Read file as Pairs object
  x = rtracklayer::import(f, format='bedpe')
  
  ## Update metadata column names
  x.mcol.names = colnames(read.csv(f, h=T, stringsAsFactors=F, sep='\t', check.names=F))
  colnames(mcols(x))[3:ncol(mcols(x))] = x.mcol.names[11:length(x.mcol.names)]
  mcols(x)$type = mcols(x)$name
  
  
  ## Convert to breakpoint ranges
  x  = StructuralVariantAnnotation::pairs2breakpointgr(x)
  
  return(x)
  
}



## Read a headered, tab-delimited CNV file into a GRanges object
readCNV = function(f) {
  
  x = read.csv(f, h=F, stringsAsFactors=F, sep='\t', comment.char='#')
  colnames(x)[1:3] = c('chr','start','end')
  x = makeGRangesFromDataFrame(x)
  
  return(x)
  
}



## Find the nearest copy number changepoints to each breakend
annotateWithClosestChangepoint = function(sv, cnv) {
  
  sv$cnv = ''
  cnv.str = paste0(as.character(seqnames(cnv)), ':', start(cnv), '-', end(cnv))
  nearest.cnv = GenomicRanges::nearest(sv, cnv)
  
  ## Make sure NAs (i.e. no nearest neightbor) are preserved as blanks  
  idx.na = which(is.na(nearest.cnv))
  nearest.cnv[idx.na] = 1  
  
  sv$cnv = cnv.str[nearest.cnv]
  sv$cnv[idx.na] = '' 
  
  return(sv)  
  
}



## Convert breakpointRanges to BEDPE
vcfToBedpe = function(vcf) {
  
  sqn = as.character(seqnames(vcf))
  strand = as.character(strand(vcf))
  res = c()
  processed = c()
  
  for (i in 1:length(vcf)) {
    bnd = names(vcf)[i]
    partner = vcf$partner[i]
    partner.idx = which(names(vcf) == partner)
    
    ## If we don't have exactly one partner, exclude this variant
    if (length(partner.idx) != 1) {
      warning('Missing partner for breakend ', bnd)
      next
    }
    
    ## Check to see if we've alrady processed this or it's partner
    if (any(c(bnd, partner) %in% processed)) {
      next
    }
    
    
    ## Combine breakends in single line
    res.i = c(sqn[i], start(vcf)[i], end(vcf)[i],                                  ## chr1, start1, end1
              sqn[partner.idx], start(vcf)[partner.idx], end(vcf)[partner.idx],    ## chr2, start2, end 2
              vcf$type[i], '.', strand[i], strand[partner.idx],                    ## type, score, strand1, strand2
              vcf$evidence[i], vcf$tools[i], vcf$`tumor--normal`[i], vcf$info[i],  ## evidence, tools, TN, info 
              vcf$cnv[i], vcf$cnv[partner.idx])                                    ## changepoint1 , changpoint2
    ## Add to result, keep track of processed breakends
    res = rbind(res, res.i)
    processed = c(processed, bnd, partner)
  } 
  
  
  ## Add colnames and fill in simple event classifications
  colnames(res) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 
                    'score', 'strand1', 'strand2', 'evidence', 'tools', 'tumor--normal',
                    'info', 'cnv_changepoint_1', 'cnv_changepoint_2')
  res = as.data.frame(res, stringsAsFactors=F)
  
  
  ## Fix coordinates (have to subtract when starting from a bedpe)
  res$start1 = as.numeric(res$start1) - 1
  res$start2 = as.numeric(res$start2) - 1
  
  
  colnames(res)[1] = paste0('#', colnames(res)[1])
  
  return(res)
  
}



## Collect arguments
option_list = list(
  make_option(c("-b", "--bedpe"),    type='character',  help="Input BEDPE"),
  make_option(c("-c", "--cnv"),      type='character',  help="BED file containing CNV intervals"),
  make_option(c("-o", "--out_file"), type='character',  help="Output BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))



## Read bedpe
sv = readBEDPE(opt$bedpe)

## Read CNV
cnv = readCNV(opt$cnv)

## Annotate breakpoints with closest changepoint
sv = annotateWithClosestChangepoint(sv=sv, cnv=cnv)

## Convert to bedpe 
res = vcfToBedpe(sv)

## Write result
write.table(res, opt$out_file, row.names=F, col.names=T, sep='\t', quote=F)
