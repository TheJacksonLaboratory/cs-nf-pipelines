suppressMessages(library(ASCAT))

### ASCAT Run ######

# Note: this script expects ASCAT is running on single sample BAF/LRR files.

args=(commandArgs(TRUE))

sampleID = args[1]
BAF_file = args[2]
LRR_file = args[3]
gender   = args[4]
platform = args[5]
GC_file  = args[6]
RT_file  = args[7]

######

# Expected SNP POS file: 
# Probe Set ID    Chromosome      Physical Position
# CN_473963       1       61736
# CN_473964       1       61808

## the above can be taken from the BAF file. The BAF file contains positions for all valid SNPs. 
SNPpos <- read.table(BAF_file, sep = "\t", header = TRUE)[ ,1:3]
colnames(SNPpos) <- c('Probe_Set_ID', 'Chromosome', 'Physical_Position')

##

if (gender == 'NA') {
  gender = 'XY'
}

ascat.bc = ascat.loadData(Tumor_LogR_file = LRR_file, Tumor_BAF_file = BAF_file, gender = gender, genomeVersion = "hg38")

ascat.bc$samples[1] <- sampleID
colnames(ascat.bc[["Tumor_LogR"]]) <- sampleID
colnames(ascat.bc[["Tumor_BAF"]]) <- sampleID

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")

ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = GC_file, replictimingfile = RT_file)

ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")

gg = ascat.predictGermlineGenotypes(ascat.bc, platform = platform)

ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = gg)

ascat.plotSegmentedData(ascat.bc)

ascat.output = ascat.runAscat(ascat.bc, write_segments = T)

##

QC = ascat.metrics(ascat.bc, ascat.output)

write.table(as.data.frame(QC), file = paste0(sampleID, "_sample.QC.txt"), sep="\t", quote=F, row.names=F, col.names=T)

save(ascat.bc, ascat.output, QC, file = paste0(sampleID, "_ASCAT_objects.Rdata"))

##

if ( length(ascat.output$failedarrays) == 0 ) {
  
  num_probes <- vector(mode="numeric", length=nrow(ascat.output$segments_raw))
  for (i in 1:nrow(ascat.output$segments_raw)) {
    L1 = which(SNPpos$Chromosome == ascat.output$segments_raw$chr[i] & SNPpos$Physical_Position == ascat.output$segments_raw$startpos[i])
    L2 = which(SNPpos$Chromosome ==  ascat.output$segments_raw$chr[i] & SNPpos$Physical_Position == ascat.output$segments_raw$endpos[i])
    num_probes[i] = L2[length(L2)] - L1[1] + 1
  }
  seg_raw = cbind(ascat.output$segments_raw,num_probes)
  
  num_probes <- vector(mode="numeric", length=nrow(ascat.output$segments))
  for (i in 1:nrow(ascat.output$segments)) {
    
    #print(i)
    L1 = which(SNPpos$Chromosome == ascat.output$segments$chr[i] & SNPpos$Physical_Position == ascat.output$segments$startpos[i])
    L2 = which(SNPpos$Chromosome ==  ascat.output$segments$chr[i] & SNPpos$Physical_Position == ascat.output$segments$endpos[i])
    num_probes[i] = L2[length(L2)] - L1[1] + 1
    
  }
  seg = cbind(ascat.output$segments,num_probes)
  
  seg_raw_dfs <- split(seg_raw, seg_raw$sample)
  seg_dfs <- split(seg, seg$sample)
  
  for (samp in names(seg_raw_dfs)){
    write.table(seg_raw_dfs[[samp]], file = paste0(samp, ".segments_raw.txt"), sep="\t", quote=F, row.names=F)
    write.table(seg_dfs[[samp]], file = paste0(samp, ".segments.txt"), sep="\t", quote=F, row.names=F)
    write.table(as.data.frame(ascat.output$aberrantcellfraction)[row.names(as.data.frame(ascat.output$aberrantcellfraction)) %in% samp,], file=paste(samp,".aberrantcellfraction.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(as.data.frame(ascat.output$ploidy)[row.names(as.data.frame(ascat.output$ploidy)) %in% samp,], file=paste(samp,".ploidy.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  }
  
} else {
  
  write.table(as.data.frame(ascat.output$failedarrays), file="ASCAT.failedarrays.txt", sep="\t", quote=F, row.names=F, col.names=F)
  
}

if ( !is.null(ascat.output$nonaberrantarrays) ) {
  
  write.table(as.data.frame(ascat.output$nonaberrantarrays), file="ASCAT.nonaberrantarrays.txt", sep="\t", quote=F, row.names=F, col.names=F)
  
}

sessionInfo()