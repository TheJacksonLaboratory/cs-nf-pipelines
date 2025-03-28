## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Annotate a merged & annotated BEDPE with gene information
libs = c('optparse', 'StructuralVariantAnnotation', 'VariantAnnotation', 'rtracklayer', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)


CLOSEST_MAX_DISTANCE = 2e4     ## For intergenic CNVs, ignore distanceToNearest() hits farther than this



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



## Read Ensembl gene info
readEnsembl = function(f) {
  
  ## Read and get rid of unnecessary columns  
  x = read.csv(f, h=F, stringsAsFactors=F, sep='\t')
  colnames(x) = c('gene_chr', 'gene_start', 'gene_end', 'strand', 'name', 'na1', 'na2', 'exons', 'exon_starts', 'exon_ends', 'na3', 'intron_starts', 'intron_ends')
  x = x[, !grepl('na[0-9]$', colnames(x))]
  
  x = GenomicRanges::makeGRangesFromDataFrame(x, 
                                              keep.extra.columns=T, 
                                              seqnames.field='gene_chr', 
                                              start.field='gene_start',
                                              end.field='gene_end')
  
  x$intron_starts = gsub(',$','',x$intron_starts)
  x$intron_ends = gsub(',$','',x$intron_ends)
  
  return(x)
  
}



## Check if breakpoints fall within genes
annotateWithDisruptions = function(sv, genes) {
  
  genes$disrupt = genes$name
  sv = gr.val(query=sv, target=genes, val='disrupt')
  sv$disrupt = gsub(' ', '', sv$disrupt, fixed=T)
  
  return(sv)
  
}



readCancerCensus = function(f) {
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  colnames(x) = c('chrom', 'start', 'end', 'name', 'locus')
  
  
  # x$name = gsub('\\|.*$', '', x$name)
  x = x[, !colnames(x) %in% 'locus']
  
  x = GenomicRanges::makeGRangesFromDataFrame(x, 
                                              keep.extra.columns=T, 
                                              seqnames.field='chrom', 
                                              start.field='start',
                                              end.field='end')
  
  return(x)
  
}



## Check if breakends fall within introns
annotateWithIntronic = function(sv=sv, genes=genes) {
  
  sv$intronic = ''
  
  ## Subset gene list to those that have introns and are already known to overlap with breakpoints 
  genes = genes[genes$intron_starts != '-' & genes$name %in% unique(unlist(strsplit(sv$disrupt, ',')))]
  
  ## Expand introns to GRanges object (keeping track of what genes they belong to)
  introns = base::mapply(function(x, y, z, n) GRanges(x, IRanges(as.numeric(y), as.numeric(z), name=rep(n, length(y)))),
                         x=as.character(seqnames(genes)),
                         y=strsplit(genes$intron_starts, ',', fixed=T),
                         z=strsplit(genes$intron_ends, ',', fixed=T),
                         n=genes$name)
  
  introns = Reduce(c, introns)
  
  
  ## Do SV breakends overlap any introns? 
  processed = c()
  for (i in 1:length(sv)) {
    partner.idx = which(names(sv) == sv$partner[i])
    
    ## Check if we've already procesed the full breakpoint
    if (i %in% processed) {
      next
    }    
    
    ## Check which introns these breakpoints overlap
    intron.hits = findOverlaps(query=c(sv[c(i, partner.idx)]), subject=introns)
    
    ## Do we have any duplicated intron hits? If so, these breakends fall within the same intron
    gene.names = names(introns)[subjectHits(intron.hits)[duplicated(subjectHits(intron.hits))]]
    sv$intronic[c(i, partner.idx)] = paste(gene.names, collapse=',')      
    
    processed = c(processed, i, partner.idx)
    
  }
  
  return(sv)
  
}



## Find the nearest copy number changepoints to each breakend
annotateWithClosest = function(sv, genes, closest.max.distance=CLOSEST_MAX_DISTANCE) {
  
  sv$closest = ''
  
  ## For each breakend
  for (i in 1:length(sv)) {
    
    ## Find the nearest non-overlapping gene, i.e.
    ## exclude any disrupted/contained gene(s) from the comparisons
    disrupt = unlist(strsplit(sv$disrupt[i], ','))
    # contains = unlist(strsplit(sv$contains[i], ','))
    intronic = unlist(strsplit(sv$intronic[i], ','))
    
    if (length(c(disrupt, intronic)) > 0) {
      genes.i = genes[-which(genes$name %in% c(disrupt, intronic))]
    } else {
      genes.i = genes
    }   
    

    ## Add closest gene, subject to distance cutoff (default 20kb)
    closest = GenomicRanges::distanceToNearest(x=sv[i], subject=genes.i, ignore.strand=T)
    closest = closest[mcols(closest)$distance <= closest.max.distance]
    
    if (length(closest) > 0) {
      sv$closest[i] = genes.i$name[subjectHits(closest)]  
    }
 
  }
  
  return(sv)  
  
}



annotateWithContained = function(sv, genes, sv.colname='contains', allow.partial.overlap=F) {
  
  mcols(sv)[,sv.colname] = ''
  
  ## For each breakend
  processed = c()
  for (i in 1:length(sv)) {
    
    partner.idx = which(names(sv) == sv$partner[i])
    is.translocation = as.character(seqnames(sv[i])) != as.character(seqnames(sv[partner.idx]))
    contains = '' 
    
    ## Check if we've already procesed the full breakpoint
    if (i %in% processed) {
      next
    }    
    
    ## Only check contained genes if this is intrachromosomal
    if (!is.translocation) {
      
      ## If allow.partial.overlap=T, allow partially overlapping intervals, otherwise, just a 
      ## simple coordinate check to see if intervals are fully contained
      if (allow.partial.overlap) {
      
        bkpt = GRanges(as.character(seqnames(sv))[i], IRanges(start(sv)[i], end(sv)[partner.idx]))  
        contains = genes$name[genes %^% bkpt]
              
      } else {
        
        contains = genes$name[start(sv)[i] <= start(genes) & 
                                end(genes) <= end(sv)[partner.idx] & 
                                as.character(seqnames(genes)) == as.character(seqnames(sv[i]))]          
      
      }
      
      
      
    } else if (allow.partial.overlap && is.translocation) {
      
      contains = genes$name[genes %^% sv[c(i, partner.idx)]]
      
    } 
    
    mcols(sv)[i, sv.colname] = paste(contains, collapse=',')
    processed = c(processed, i, partner.idx)
    
  }
  
  return(sv)
  
}



## Convert breakpointRanges to BEDPE
vcfToBedpe = function(vcf, supplemental) {
  
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
    
    
    ## Add disrupt/closest/intronic to info field
    disrupt.l = paste0('DisruptL=', vcf$disrupt[i])
    disrupt.r = paste0('DisruptR=', vcf$disrupt[partner.idx])
    closest.l = paste0('ClosestL=', vcf$closest[i])
    closest.r = paste0('ClosestR=', vcf$closest[partner.idx])
    intronic  = paste0('Intronic=', vcf$intronic[i])
    contained = paste0('Contained=', vcf$contains[i])
    
    if (supplemental) {
      gene.str = paste(disrupt.l, disrupt.r, closest.l, closest.r, intronic, contained, sep=';')
    } else {
      gene.str = paste(disrupt.l, disrupt.r, closest.l, closest.r, intronic, sep=';')  
    }
    
    
    
    ## Add cancer census genes only if they exist
    if (vcf$cgc[i] != '' || vcf$cgc[partner.idx] !='') {
      cgc.str = paste0('Cancer_census=', paste(setdiff(vcf$cgc[c(i,partner.idx)],''),collapse=','))
      gene.str = paste(gene.str, cgc.str, sep=';')
    }
    
    vcf$info[i] = paste0(vcf$info[i], gene.str)
    
    
    ## Combine breakends in single line
    res.i = c(sqn[i], start(vcf)[i], end(vcf)[i],                                  ## chr1, start1, end1
              sqn[partner.idx], start(vcf)[partner.idx], end(vcf)[partner.idx],    ## chr2, start2, end 2
              vcf$type[i], '.', strand[i], strand[partner.idx],                    ## type, score, strand1, strand2
              vcf$evidence[i], vcf$tools[i], vcf$`tumor--normal`[i], vcf$info[i])  ## evidence, tools, TN, info 
    
    ## Add to result, keep track of processed breakends
    res = rbind(res, res.i)
    processed = c(processed, bnd, partner)
  } 
  
  
  ## Add colnames and fill in simple event classifications
  colnames(res) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 
                    'score', 'strand1', 'strand2', 'evidence', 'tools', 'tumor--normal',
                    'info')
  res = as.data.frame(res, stringsAsFactors=F)
  
  ## Fix coordinates (have to subtract when starting from a bedpe)
  res$start1 = as.numeric(res$start1) - 1
  res$start2 = as.numeric(res$start2) - 1
  
  
  colnames(res)[1] = paste0('#', colnames(res)[1])
  
  return(res)
  
}



## Collect arguments
option_list = list(
  make_option(c("-b", "--bedpe"),         type='character',    help="Input BEDPE"),
  make_option(c("-e", "--ensembl"),       type='character',    help="Ensembl gene list"),
  make_option(c("-c", "--cancer_census"), type='character',    help="Cancer census gene list"),
  make_option(c("-s", "--supplemental"),  action='store_true', default=F,  help="Add supplementary gene annotations?"),
  make_option(c("-o", "--out_file"),      type='character',    help="Output BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))


## Read bedpe
sv <- tryCatch( 
  {
      readBEDPE(opt$bedpe)
  },
  error = function(e) {
      res <- data.frame('a'=character(), 'b'=numeric(), 'c'=numeric(), 'd'=character(), 'e'=numeric(), 'f'=numeric(), 'g'=character(), 'h'=character(), 'i'=character(), 'j'=character(), 'k'=character(), 'l'=character(), 'm'=character(), 'n'=character())
      colnames(res) = c('#chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'score', 'strand1', 'strand2', 'evidence', 'tools', 'tumor--normal', 'info')
      write.table(res, opt$out_file, row.names=F, col.names=T, sep='\t', quote=F)
      quit(save = "no", status = 0)
  }
)
# if SV CSV is empty (i.e., no somatic SV), write an empty file as output object.

## Read gene lists
genes = readEnsembl(opt$ensembl)
cgc = readCancerCensus(opt$cancer_census)

## Add contained ensembl and cgc genes
sv = annotateWithContained(sv=sv, genes=genes, allow.partial.overlap=F)
sv = annotateWithContained(sv=sv, genes=cgc, sv.colname='cgc', allow.partial.overlap=T)

## Add disrupted genes
sv = annotateWithDisruptions(sv=sv, genes=genes)

## Which breakpoints fall within introns?
sv = annotateWithIntronic(sv=sv, genes=genes)

## Add closest (non-contained) genes
sv = annotateWithClosest(sv=sv, genes=genes)

## Convert to bedpe 
res = vcfToBedpe(sv, supplemental=opt$supplemental)

## Write result
write.table(res, opt$out_file, row.names=F, col.names=T, sep='\t', quote=F)
