## Filter a bedpe for somatic variants (i.e., not in specified germline databases), and 
## high-confidence variants (2+ callers or 1 caller with a nearby copy number changepoint)
libs = c('optparse', 'GenomicRanges')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)


## Check if databases db are in info string x
inDatabase = function(x, db) {
  
  ## Split info field, look for database entry
  x = unlist(strsplit(x, ';', fixed=T))
  x = grep('known=', x, fixed=T, value=T)
  x = gsub('known=', '', x, fixed=T)
  x = unlist(strsplit(x, ',', fixed=T))
  
  return(any(x %in% db))
  
}



makeGRangesFromChangepoint = function(x) {
  
  x = unlist(strsplit(x, ':|-'))
  GRanges(seqnames=x[1], ranges=IRanges(as.numeric(x[2:3]), as.numeric(x[2:3])))
                                              
}



## Is variant x a high-confidence variant? 
## Meant to be used with apply(,2,)
isHighConfidence = function(x, cpmax) {
  
  ## Is there support from multiple callers?
  multi.caller = grepl(',', x['tools'])

  ## Is either breakpoint close enough to its nearest changepoint?
  x1.gr = GRanges(seqnames=x['#chr1'], ranges=IRanges(as.numeric(x['start1']), as.numeric(x['end1'])))
  ch1.gr = makeGRangesFromChangepoint(x['cnv_changepoint_1'])
  near.ch1 = any(GenomicRanges::distance(x1.gr, ch1.gr) <= cpmax)
  
  x2.gr = GRanges(seqnames=x['chr2'], ranges=IRanges(as.numeric(x['start2']), as.numeric(x['end2'])))
  ch2.gr = makeGRangesFromChangepoint(x['cnv_changepoint_2'])
  near.ch2 = any(GenomicRanges::distance(x2.gr, ch2.gr) <= cpmax)
  
  return(multi.caller || near.ch1 || near.ch2)
  
}



## Collect arguments
option_list = list(
  make_option(c("-b", "--bedpe"),                     type='character',  help="Input BEDPE"),
  make_option(c("-m", "--max_changepoint_distance"),  type='numeric',    help="Maximum distance a changepoint can be from a breakpoint to 'rescue' it into the high-confidence set"),
  make_option(c("-f", "--filter_databases"),          type='character',  help="Comma-separated list of databases to filter, looking in the info field"),
  make_option(c("-s", "--out_file_somatic"),          type='character',  help="Output somatic BEDPE"),
  make_option(c("-o", "--out_file_highconf"),         type='character',  help="Output high-confidence BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))

## Unpack arguments
opt$filter_databases = unlist(strsplit(opt$filter_databases, ',', fixed=T))


## Read bedpe, filter for known germline variants 
x = read.csv(opt$bedpe, h=T, stringsAsFactors=F, sep='\t', check.names=F)
x = x[!sapply(x$info, inDatabase, opt$filter_databases), ]

## Write out somatic variants
write.table(x, opt$out_file_somatic, row.names=F, col.names=T, sep='\t', quote=F)

## Filter for high confidence
x = x[apply(x, 1, isHighConfidence, opt$max_changepoint_distance), ]

## Write result
write.table(x, opt$out_file_highconf, row.names=F, col.names=T, sep='\t', quote=F)
