## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Merge arbitrary number of VCFs, annotate with simple event type 
libs = c('optparse', 'StructuralVariantAnnotation', 'VariantAnnotation', 'rtracklayer', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)


SUPPORTED_CALLERS = c('manta', 'lumpy', 'svaba', 'gridss', 'delly')     ## Update this flag when adding support for new callers
SVABA_MIN_LENGTH = 1001                                        ## Svaba-unique calls shorter than this appear to be artifactual


## Callers have different names for the same pieces of evidence,
## For now handle each case separately
## TODO: Add support for GRIDSS 
getReadSupport = function(vcf, caller, sample_id, supplementary=FALSE, supported_callers=SUPPORTED_CALLERS) {
  
  ## Don't try to process genotype info if we don't know how
  if (!caller %in% supported_callers) {
    stop('Caller ', caller, ' is not currently supported. Supported callers: ', paste(supported_callers, collapse=','))
  }
  
  ## It's a possibility that the sample names in the VCF will be 
  ## the full path to the BAM used instead of just the sample ID
  ## Just grab the index of the correct column
  if (!sample_id %in% colnames(geno(vcf)[[1]])) {
    sample_id = which(gsub('\\.final\\.bam$','',basename(colnames(geno(vcf)[[1]]))) %in% sample_id)
  }
  
  
  if (caller == 'manta') {
    
    ## Common info
    sr = geno(vcf)$SR[, sample_id]
    sr = sapply(sr, `[`, 2)
    pe = geno(vcf)$PR[, sample_id]
    pe = sapply(pe, `[`, 2)
    
    ## Supplementary info
    supp_string = paste0(caller,'_SOMATICSCORE=', info(vcf)$SOMATICSCORE) 
    
  } else if (caller == 'svaba') {
    
    ## Common info
    sr = geno(vcf)$SR[, sample_id]
    pe = geno(vcf)$DR[, sample_id]
    
    ## Supplementary info
    ad = paste0(caller,'_AD=', geno(vcf)$AD[, sample_id])
    dp = paste0(caller,'_DP=', geno(vcf)$DP[, sample_id])
    lo = paste0(caller,'_LO=', geno(vcf)$LO[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    supp_string = paste(ad, dp, lo, gt, sep=',')
    
  } else if (caller == 'lumpy') {
    
    ## Common info
    sr = info(vcf)$SR
    pe = info(vcf)$PE
    ## Supplementary info
    ro = paste0(caller,'_RO=', geno(vcf)$RO[, sample_id])
    ao = paste0(caller,'_AO=', geno(vcf)$AO[, sample_id])
    dp = paste0(caller,'_DP=', geno(vcf)$DP[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    supp_string = paste(ro, ao, dp, gt, sep=',')
    
    
  } else if (caller == 'gridss') {
    
    ## Common info
    sr = geno(vcf)$SR[, sample_id]
    pe = geno(vcf)$RP[, sample_id]
    
    ## Supplementary info
    vf = paste0(caller,'_VF=', geno(vcf)$VF[, sample_id])
    asq = paste0(caller,'_ASQ=', geno(vcf)$ASQ[, sample_id])
    qual = paste0(caller,'_QUAL=', geno(vcf)$QUAL[, sample_id])
    supp_string = paste(vf, asq, qual, sep=',')
    
  } else if (caller == 'delly') {
    ## Common info
    sr = info(vcf)$SR
    # sr = sapply(sr, `[`, 2)
    pe = info(vcf)$PE
    # pe = sapply(pe, `[`, 2)
    ## Supplementary info
    dr = paste0(caller,'_DR=', geno(vcf)$DR[, sample_id])
    dv = paste0(caller,'_DV=', geno(vcf)$DV[, sample_id])
    rr = paste0(caller,'_RR=', geno(vcf)$RR[, sample_id])
    rv = paste0(caller,'_RV=', geno(vcf)$RV[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    supp_string = paste(dr, dv, rr, rv, gt, sep=',')
  }
  
  ## Set NA to 0
  ## TODO: Keep this? 
  sr[is.na(sr)] = 0
  pe[is.na(pe)] = 0
  
  ## Build output string
  if (supplementary) {
    res = paste0('[',caller,'_SR=',sr,',', caller,'_PE=', pe,',', supp_string,']')
  } else {
    res = paste0('[',caller,'_SR=',sr,',', caller,'_PE=', pe,']')  
  }
  
  return(res)
  
}



sumSupport = function(x) {
  sapply(str_extract_all(x, '(?<=\\=)[0-9]+(?=,|\\])'), function(y) sum(as.numeric(y)))
}



removeRedundantBreakpoints = function(x) {
  
  ## Find duplicates
  key = unlist(strsplit(x$breakendPosID,','))
  key.count = table(key)
  key.dup = key.count[key.count > 1]
  
  
  ## If there aren't duplicates we don't have anything to do
  if (length(key.dup) == 0) {
    return(x)
  }
  
  
  ## For each set of duplicate breakends, select the one with the higher score
  x.idx.rm = c()
  for (i in names(key.dup)) {
    
    ## Subset to breakends of interest
    x.idx = grep(i, x$breakendPosID, fixed=T)
    xi = x[x.idx]
    
    ## Collect support
    xi$read.support = sumSupport(xi$support)
    xi$multicaller.support = grepl('],[',xi$support,fixed=T)
    
    
    
    ## Automatically keep breakends with multi-caller support
    idx.multi = which(xi$multicaller.support)
    if (length(idx.multi) > 0) {
      x.idx.rm = c(x.idx.rm, x.idx[-idx.multi])
      next
    }
    
    
    ## Automatically discard breakends with the lowest support
    idx.max = which(xi$read.support %in% max(xi$read.support))
    if (length(idx.max) > 0) {
      x.idx.rm = c(x.idx.rm, x.idx[-idx.max])     
      x.idx = x.idx[idx.max]
      xi = xi[idx.max]
    }
    
    
    ## If there are multiple breakends tied for highest read support
    if (length(xi) > 1) {
      
      if (all(!is.na(xi$svLen))) {
        
        ## If all are non-TRA take longest SV 
        x.idx.rm = c(x.idx.rm, x.idx[-which.max(abs(xi$svLen))])
        
      } else if (all(is.na(xi$svLen)) && length(unique(as.character(seqnames(xi)))) == 1) {
        
        ## If all TRA to the same chr select rightmost coordinate
        partners = x[names(x) %in% xi$partner]
        partner.keep = names(partners)[which.max(start(partners))]
        x.idx.rm = c(x.idx.rm, x.idx[!xi$partner %in% partner.keep])
        
      }
      
      ## Otherwise, just keep tied SVs
      
    }
    
  }
  
  
  ## Remove breakends and their partners if we have any to remove
  if (length(x.idx.rm) > 0) {
    x = x[-x.idx.rm]  
    x = x[names(x) %in% x$partner]
  }
  
  return(x)
  
}



## Compute error between query and subject for a hits object
computeError = function(query, subject, hits) {
  
  ## Init result dataframe
  error = data.frame(local=rep(NA, length(queryHits(hits))), 
                     remote=rep(NA, length(queryHits(hits))))
  
  
  ## For each hit
  for (i in 1:length(queryHits(hits))) {
    
    ## Compute local error (error between breakends at hit i) and remote error (error between the partners of 
    ## the breakends at hit i)
    error$local[i] = StructuralVariantAnnotation:::.distance(query[queryHits(hits)[i]], subject[subjectHits(hits)[i]])$min
    error$remote[i] = StructuralVariantAnnotation:::.distance(query[names(query) == query[queryHits(hits)[i]]$partner],
                                                              subject[names(subject) == subject[subjectHits(hits)[i]]$partner]
                                                              )$min
  }
  
  return(error)
    
}



## Take the union of callsets a and b, both breakpointRanges objects
## If multiple hits found in b for a, choose the closest match, measured
## as the mean distance between breakends
mergeCallsets = function(a, b, slop) {
  
  ## Find overlaps
  overlaps = StructuralVariantAnnotation::findBreakpointOverlaps(query=a, 
                                                                subject=b, 
                                                                maxgap=slop, 
                                                                sizemargin=0.8, 
                                                                restrictMarginToSizeMultiple=0.8)
  
  
  ## If we have any duplicate query hits, choose hit based on match quality
  if(anyDuplicated(queryHits(overlaps))) {
    
    ## Compute local and remote breakend basepair error on matches  
    error = computeError(query=a, subject=b, hits=overlaps)
    ## Get duplicate hits  
    dup.query.hits = table(queryHits(overlaps))
    dup.query.hits = names(dup.query.hits[dup.query.hits > 1])
    
    ## Determine which hits we're removing
    idx.hits.rm = c()
    for (d in dup.query.hits) {
      
      idx.dup.query.hits = which(queryHits(overlaps) %in% d)
      local.error = error$local[idx.dup.query.hits]
      remote.error = error$remote[idx.dup.query.hits]
      mean.error = rowMeans(cbind(local.error, remote.error))
      
      ## Keep the query hit with the smallest mean error
      idx.hits.rm = c(idx.hits.rm, idx.dup.query.hits[which.max(mean.error)])
      
    }
    
    overlaps = overlaps[-idx.hits.rm]   
  }
  
  ## For matching SVs, merge caller support
  a$support[queryHits(overlaps)] = paste0(a$support[queryHits(overlaps)],',',b$support[subjectHits(overlaps)])
  a$breakendPosID[queryHits(overlaps)] = paste0(a$breakendPosID[queryHits(overlaps)],',',b$breakendPosID[subjectHits(overlaps)])
  
  ## Pull in non-matching SVs from b
  res = c(a, b[-subjectHits(overlaps)])

  return(res)
  
}



## Convert breakpointRanges to BEDPE
vcfToBedpe = function(vcf, supplemental=F) {
  
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
    
    ## Which support column should we use? 
    if (supplemental) {
      support = vcf$supplemental[i] 
    } else {
      support = vcf$support[i] 
    }
    
    
    ## Combine breakends in single line
    res.i = c(sqn[i], start(vcf)[i], end(vcf)[i],                                  ## chr1, start1, end1
              sqn[partner.idx], start(vcf)[partner.idx], end(vcf)[partner.idx],    ## chr2, start2, end 2
              'BND', '.', strand[i], strand[partner.idx], support)                 ## type, score, strand1, strand2, support
    
    ## Add to result, keep track of processed breakends
    res = rbind(res, res.i)
    processed = c(processed, bnd, partner)
  } 
  
  
  ## Add colnames and fill in simple event classifications
  colnames(res) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'score', 'strand1', 'strand2', 'evidence')
  res = as.data.frame(res, stringsAsFactors=F)
  
  res$type[res$strand1 == '+' & res$strand2 == '-'] = 'DEL'
  res$type[res$strand1 == '-' & res$strand2 == '+'] = 'DUP'
  res$type[res$strand1 == '-' & res$strand2 == '-'] = 'INV'
  res$type[res$strand1 == '+' & res$strand2 == '+'] = 'INV'
  res$type[res$chr1 != res$chr2] = 'TRA'
  
  ## Sort by chromosome 
  res = res[order(factor(res$chr1, levels=levels(seqnames(vcf))), res$start1, res$end1, decreasing=F), ]
  
  ## Simplify coordinates
  res$end1 = as.numeric(res$start1) + 1
  res$end2 = as.numeric(res$start2) + 1
  
  
  ## Extract tool info from read support column
  res$tools = sapply(res$evidence, function(x) paste(unlist(stringr::str_extract_all(x, '(?<=\\[)[a-z]+(?=_)')), collapse=','))
  
  colnames(res)[1] = paste0('#', colnames(res)[1])
  
  return(res)
  
}



## Collect arguments
option_list = list(
  make_option(c("-v", "--vcf"),                   type='character', help="Comma-delimited list of breakend notation VCFs"),
  make_option(c("-c", "--callers"),               type='character', help="Comma-delimited list of SV caller names corresponding to the order of VCFs given in --vcf"),
  make_option(c("-t", "--tumor"),                 type='character', help="Tumor sample ID"),
  make_option(c("-n", "--normal"),                type='character', help="Normal sample ID"),
  make_option(c("-b", "--build"),                 type='character', help="Genome build"),
  make_option(c("-s", "--slop"),                  type='numeric',   help="Padding to use when comparing breakpoints"),
  make_option(c("-l", "--min_sv_length"),         type='numeric',   help="Filter SVs shorter than this length"),
  make_option(c("-a", "--allowed_chr"),           type='character', help="Comma-delimited list of chromosomes to keep"),
  make_option(c("-o", "--out_file"),              type='character', help="Output BEDPE"), 
  make_option(c("-p", "--out_file_supplemental"), type='character', help="Output supplemental BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))



## Unpack arguments
opt$vcf = unlist(strsplit(opt$vcf, ',', fixed=T))
opt$callers = unlist(strsplit(opt$callers, ',', fixed=T))
opt$allowed_chr = unlist(strsplit(opt$allowed_chr, ',', fixed=T))



## Iteratively merge VCFs
res = NULL
for (i in 1:length(opt$vcf)) {
  ## Read VCF
  caller = opt$caller[i]
  vcf = VariantAnnotation::readVcf(opt$vcf[i], genome=opt$build)
  
  ## Get read support
  rowRanges(vcf)$support = getReadSupport(vcf=vcf, caller=caller, sample_id=opt$tumor)
  rowRanges(vcf)$supplemental = getReadSupport(vcf=vcf, caller=caller, sample_id=opt$tumor, supplementary=T )
  
  ## Convert to breakpointRanges object, don't adjust for CIPOS uncertainty (i.e. keep nominalPosition). 
  ## For Manta, infer missing breakpoint is required as the caller does not insert the recip call in the VCF as the other calls do. 
  vcf = StructuralVariantAnnotation::breakpointRanges(vcf, nominalPosition=T, inferMissingBreakends=T)
  ## Add breakendPosID for later redundancy checks
  vcf$breakendPosID = paste0('[',caller,'=',as.character(seqnames(vcf)),':',start(vcf),':',strand(vcf),']')
  ## Overlap if this isn't the first callset
  if (i == 1) {
    res = vcf
  } else {
    res = mergeCallsets(a=res, b=vcf, slop=opt$slop)
  }
}

## Handle breakpoints with duplicate start or end positions
res = removeRedundantBreakpoints(res)

## Convert to bedpe, apply some filters 
for (i in c('main','supplemental')) {
  
  outfile = ifelse(i=='main', opt$out_file, opt$out_file_supplemental)
  
  ## Convert to BEDPE format
  res.i = vcfToBedpe(res, supplemental=i=='supplemental')
  res.i$`tumor--normal` = paste0(opt$tumor,'--',opt$normal)
  
  ## Filter non-TRA variants for minimum length opt$min_sv_length
  sv.lengths = abs(as.numeric(res.i$start2) - as.numeric(res.i$start1))
  res.i = res.i[res.i$type == 'TRA' | sv.lengths >= opt$min_sv_length, ]
  
  ## Filter non-TRA svaba-unique variants less than SVABA_MIN_LENGTH
  sv.lengths = abs(as.numeric(res.i$start2) - as.numeric(res.i$start1))
  res.i = res.i[(res.i$tools != 'svaba' | res.i$type == 'TRA') | (res.i$tools == 'svaba' & sv.lengths >= SVABA_MIN_LENGTH), ]
  
  ## Filter SVs not occurring in allowed chromosomes (i.e. autosomes and sex chromosomes)
  res.i = res.i[res.i$`#chr1` %in% opt$allowed_chr & res.i$chr2 %in% opt$allowed_chr, ]
  
  ## Write result
  write.table(res.i, outfile, row.names=F, col.names=T, sep='\t', quote=F)
  
}
