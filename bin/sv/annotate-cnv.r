## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Annotate a merged bedpe with arbitrary databases
libs = c('optparse', 'gUtils', 'GenomicRanges', 'rtracklayer')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)

## TODO: Move to config? 
CLOSEST_MAX_DISTANCE = 2e4     ## For intergenic CNVs, ignore nearest() hits farther than this
LARGESCALE_MIN = 3e6           ## Any events smaller than this are considered focal 
DUP_LOG2 = 0.2                 ## log2 ratio cutoff for considering an event a duplication
DEL_LOG2 = -0.235              ## log2 ratio cutoff for considering an event a deletion

## Read BIC-Seq2 output into a GRanges object
## Optionally subset to CNVs in chr
readCNV = function(f, chr=NULL) {
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  colnames(x)[colnames(x) == 'log2.copyRatio'] = 'log2'
  
  x = GenomicRanges::makeGRangesFromDataFrame(x, 
                                              keep.extra.columns=T, 
                                              seqnames.field='chrom', 
                                              start.field='start',
                                              end.field='end')
  
  
  if (!is.null(chr)) {
    x = x[as.character(seqnames(x)) %in% chr]  
  }
  
  return(x)
  
}

## Read cytoband into a GRanges object
readCytoband = function(f) {
  
  x = read.csv(f, h=F, stringsAsFactors=F, sep='\t')
  colnames(x) = c('chrom', 'start', 'end', 'cytoband', 'stain')
  
  x = x[, !colnames(x) %in% 'stain']
  
  x = GenomicRanges::makeGRangesFromDataFrame(x, 
                                              keep.extra.columns=T, 
                                              seqnames.field='chrom', 
                                              start.field='start',
                                              end.field='end')
  
  return(x)
  
}

readDB = function(f) {

  x <- import(f, format = 'BED')

  return(x)
  
}

readCancerCensus = function(f) {
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  colnames(x) = c('chrom', 'start', 'end', 'cgc', 'locus')
  
  # x$cgc = gsub('\\|.*$', '', x$cgc)
  x = x[, !colnames(x) %in% 'locus']
  
  x = GenomicRanges::makeGRangesFromDataFrame(x, 
                                              keep.extra.columns=T, 
                                              seqnames.field='chrom', 
                                              start.field='start',
                                              end.field='end')
  
  return(x)
  
}

## Read Ensembl 
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

## Simplify comma-delimited cytoband list to only the first and last cytobands
.simplifyCytoband = function(x, delim=', ', collapse='-') {
  
  x = unlist(strsplit(x, delim, fixed=T))
  
  if (length(x) > 1) {
    x = paste0(x[1], collapse, x[length(x)])
  }
  
  return(x)
  
}

## Annotate with cytoband
annotateCytoband = function(cnv, cytoband) {
  
  ## Pull in cytoband info 
  cnv = cnv %$% cytoband
  
  ## Simplify comma-delimited representation to hyphenated if necessary
  cnv$cytoband = sapply(cnv$cytoband, .simplifyCytoband)
  
  ## Add chromosome information
  cnv$cytoband = paste0(as.character(seqnames(cnv)), cnv$cytoband)
  
  return(cnv)
  
}

## Annotate with databases, subject to reciprocal overlap criteria
annotateDB = function(x, db, name, overlap) {
  
  ## Find hits
  hits = GenomicRanges::findOverlaps(query=x, subject=db)
  
  ## Compute overlap
  mcols(hits)$intersection = width(pintersect(x[queryHits(hits)], db[subjectHits(hits)]))
  mcols(hits)$overlap_query =  mcols(hits)$intersection / width(x[queryHits(hits)])
  mcols(hits)$overlap_subject = mcols(hits)$intersection / width(db[subjectHits(hits)])
  
  ## Hits should meet the minimum reciprocal overlap cutoff
  ## To match bedtools::intersect's implementation of reciprocal oerlap, 
  ## The fraction overlap should be at least the same in each direction 
  hits = hits[mcols(hits)$overlap_query >= overlap & mcols(hits)$overlap_subject >= overlap]

  ## Annotate any hits we get 
  x$db[queryHits(hits)] = paste0(x$db[queryHits(hits)], ',', name)
  x$db = gsub('^,', '', x$db)
  
  return(x)
  
}

## Compare GRanges x to GRanges gene mcols intron_start, intron_end
.isIntronic = function(x, gene) {
  
  if (gene$intron_starts == '-') {
    
    is.intronic = F
    
  } else {
    
    introns = GRanges(as.character(seqnames(gene)), 
                      IRanges(as.numeric(unlist(strsplit(gene$intron_starts, ',', fixed=T))), 
                              as.numeric(unlist(strsplit(gene$intron_ends, ',', fixed=T)))))
    
    is.intronic = any(start(introns) <= start(x) && end(x) <= end(introns))
    
  }
  
  return(is.intronic)
  
}

## Annotate with ensembl genes 
annotateEnsembl = function(x, ens, closest.max.distance=CLOSEST_MAX_DISTANCE) {
  
  ## Init empty columns
  mcols(x)[, c('disrupt.l', 'disrupt.r', 'contains', 'intronic', 'intergenic', 'closest')] = ''
  
  ## Find hits 
  hits = GenomicRanges::findOverlaps(query=x, subject=ens)

  ## Check contains, disruption on each side, intronic  
  x.start = start(x[queryHits(hits)])
  x.end   = end(x[queryHits(hits)])
  gene.start = start(ens[subjectHits(hits)])
  gene.end = end(ens[subjectHits(hits)])
  
  contains = x.start <= gene.start & gene.end <= x.end    
  disrupt.r = gene.start <= x.end & x.end <= gene.end
  disrupt.l = gene.start <= x.start & x.start <= gene.end
  intronic = F
    
  
  ## We only need to check CNVs in introns if they don't contain their hit
  ## and intersect the gene on both CNV ends
  for (i in which(!contains & disrupt.r & disrupt.l)) {

    message('Checking potential intronic variant...')
    intronic[i] = .isIntronic(x=x[queryHits(hits[i])], gene=ens[subjectHits(hits[i])])

  }
    
  ## Concatenate genes and store
  ## The tapply() aggregates genes by query index (i.e. x index), which we use to map back to x
  contains = tapply(ens[subjectHits(hits)]$name[contains], queryHits(hits)[contains], paste, collapse=',')
  x$contains[as.numeric(names(contains))] = contains
  
  disrupt.r = tapply(ens[subjectHits(hits)]$name[disrupt.r], queryHits(hits)[disrupt.r], paste, collapse=',')
  x$disrupt.r[as.numeric(names(disrupt.r))] = disrupt.r
  
  disrupt.l = tapply(ens[subjectHits(hits)]$name[disrupt.l], queryHits(hits)[disrupt.l], paste, collapse=',')
  x$disrupt.l[as.numeric(names(disrupt.l))] = disrupt.l
  
  intronic = tapply(ens[subjectHits(hits)]$name[intronic], queryHits(hits)[intronic], paste, collapse=',')
  x$intronic[as.numeric(names(intronic))] = intronic
  
  
  ## Add intergenic and closest gene 
  x$intergenic[setdiff(1:length(x), queryHits(hits))] = 'yes'
  
  
  ## Add closest gene, subject to distance cutoff
  dist.to.nearest = GenomicRanges::distanceToNearest(x=x, subject=ens)
  dist.to.nearest = dist.to.nearest[mcols(dist.to.nearest)$distance <= closest.max.distance]
  x$closest[queryHits(dist.to.nearest)] = ens$name[subjectHits(dist.to.nearest)]
  
  return(x)
  
}

## Collect arguments
option_list = list(
  make_option(c("-c", "--cnv"),                   type='character', help="Input CNV calles"),
  make_option(c("-a", "--caller"),                type='character', help="Name of tool used to call CNVs in --cnv (only bicseq2 is currently supported)"),
  make_option(c("-t", "--tumor"),                 type='character', help="Comma-delimited list of database names corresponding to the order in --db_files"),
  make_option(c("-n", "--normal"),                type='character', help="Comma-delimited list of database files corresponding to the order in --db_names"),
  make_option(c("-b", "--cytoband"),              type='character', help="Cytoband file: headerless tab-delimited files with chr, start, end, cytoband, stain"),
  make_option(c("-d", "--db_names"),              type='character', help="Comma-delimited list of database names corresponding to the order in --db_files"),
  make_option(c("-s", "--db_files"),              type='character', help="Comma-delimited list of database files corresponding to the order in --db_names"),
  make_option(c("-e", "--ensembl"),               type='character', help="Ensembl gene list"),
  make_option(c("-l", "--allowed_chr"),           type='character', help="Comma-delimited list of chromosomes to keep"),
  make_option(c("-g", "--cancer_census"),         type='character', help="Cancer census gene list"),
  make_option(c("-f", "--overlap_fraction"),      type='numeric',   help="Fraction that database hits must overlap query interval"),
  make_option(c("-o", "--out_file_main"),         type='character', help="Main output BED"),
  make_option(c("-p", "--out_file_supplemental"), type='character', help="Supplemental output BED"))
opt = parse_args(OptionParser(option_list=option_list))


## Unpack arguments
opt$db_names = unlist(strsplit(opt$db_names, ',', fixed=T))
opt$db_files = unlist(strsplit(opt$db_files, ',', fixed=T))
opt$allowed_chr = unlist(strsplit(opt$allowed_chr, ',', fixed=T))


## Read files
cnv = readCNV(opt$cnv, chr=opt$allowed_chr)
cyto = readCytoband(opt$cytoband)
cgc = readCancerCensus(opt$cancer_census)
ensembl = readEnsembl(opt$ensembl)

## Add cytoband annotation
cnv = annotateCytoband(cnv=cnv, cytoband=cyto)

## Add tumor-normal id, caller info
cnv$`tumor--normal` = paste0(opt$tumor,'--',opt$normal)
cnv$tool = opt$caller


## Annotate focal/large-scale
cnv$focal = ifelse(width(cnv) < LARGESCALE_MIN, 'yes', 'no')


## Annotate dup/del/neu
cnv$type = 'NEU'
cnv$type[cnv$log2 > DUP_LOG2] = 'DUP'
cnv$type[cnv$log2 < DEL_LOG2] = 'DEL'


## Annotate with databases
cnv$db = ''
for (i in 1:length(opt$db_names)) {

  db.name = opt$db_names[i]
  db.file = opt$db_files[i]

  print(db.name)

  db = readDB(db.file)
  cnv = annotateDB(x=cnv, db=db, name=db.name, overlap=opt$overlap_fraction)

}

## Annotate with CGC genes
cnv = cnv %$% cgc
cnv$cgc = gsub(' ', '', cnv$cgc)

## Annotate with Ensembl genes
cnv = annotateEnsembl(x=cnv, ens=ensembl)

## Subtract 1 from the output start to adhere to BED standard 
start(cnv) = start(cnv) - 1

## Rename chr, convert to data frame
cnv = as.data.frame(cnv)
cnv$`#chr` = cnv$seqnames

## Build info field 
cnv$info = paste0('known=',cnv$db, ';Cancer_census=',cnv$cgc, ';DisruptL=',cnv$disrupt.l, ';DisruptR=', cnv$disrupt.r)
cnv$info[cnv$intergenic == 'yes'] = paste0(cnv$info[cnv$intergenic == 'yes'], ';Intergenic')
cnv$info[cnv$intergenic == 'yes'] = paste0(cnv$info[cnv$intergenic == 'yes'], ';Closest=', cnv$closest[cnv$intergenic == 'yes'])

## Fields included in main/supplemental are slightly different 
for (i in c('main', 'supplemental')) { 
  
  cnv.i = cnv[, c('#chr', 'start', 'end', 'type', 'log2', 'tool', 'tumor..normal', 'info', 'focal', 'cytoband')]
  colnames(cnv.i) = gsub('..', '--', colnames(cnv.i), fixed=T)
  outfile = ifelse(i == 'main', opt$out_file_main, opt$out_file_supplemental)
  
  if (i=='supplemental') {
    
    cnv.i$info = paste0(cnv$info,';Contained=',cnv$contains)
    
  }
  
  write.table(cnv.i, outfile, row.names=F, col.names=T, sep='\t', quote=F)
  
}
