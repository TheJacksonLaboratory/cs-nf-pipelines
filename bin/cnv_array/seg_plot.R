# updated to have Graces Sept '19 tweaks 

options(scipen = 999)

args=(commandArgs(TRUE))

# filename is name of the *segments_raw.extend.txt file
filename <- args[1]

sampleID <- args[2]

# outdir is the dir where png result will be written ... use "./" for current dir
outdir <- args[3]


CNS <-read.table(filename,header=T,sep="\t")

gender <- 'female'
sex <- 'female'

if (sex == "female") {
  CNS=CNS[CNS$chr!="Y",]
}

#title of plot
ploidy=round(CNS$ploidy[1], digits=2)
plottitle=paste( gsub("_","  ", sampleID), "   ploidy=",ploidy,sep="")

chromo <- unique(CNS$chr)
chromo
xx=0
y = c()
start=CNS$startpos
end=CNS$endpos

for (x in chromo) {
  start[CNS$chr == x]=start[CNS$chr == x]+xx
  end[CNS$chr == x]=end[CNS$chr == x]+xx
  tmp = CNS$endpos[CNS$chr == x]
  xx=tail(tmp,1)+xx
  y <- c(y, xx)
}

png(paste(outdir,sampleID,"_segmentsgenomeplot.copydiffploidy.png",sep=""), width=1300,height=600)

par(mar = c(5, 5, 8, 4))

val=CNS$copydiff_ploidy

# in the following stmt, family="serif" changes font to times-roman; cex.main=1.8 scales up the title font size 
# formerly also used: ylim=c(-7,20), 
plot( c(start,end), c(val,val), col="white", main=plottitle, xlab="Chromosome", ylab="Delta from Ploidy", 
      ylim=c(-8, max( c(val,val) ) ))

for (i in 1:length(start)) {
  if (CNS$LOH[i]==1) {
    polygon(c(start[i],end[i],end[i],start[i]),c(min(-7),min(-7),max(-6),max(-6)),col="lightsteelblue",border="lightsteelblue",lwd=2)
  }
  
}

segments(start,val,end,val,col="tomato",lwd=5)
abline(v=y,col="grey")
posy=c(-8)
i=1
l=0
for (x in chromo) {
  posx=(l+y[i])/2
  text(posx,posy,x,cex=1.2,srt=45)
  l=y[i]
  i=i+1
}
abline(h=0,col="black",lty=2,lwd=1.5)

# formerly also used: inset=c(0,-0.1), 
legend("topright", inset=-0.1, c("Difference from sample ploidy   ", "LOH"), xpd=TRUE, horiz=T, 
       bty="n", lty=c(1,1), lwd=6, col=c("tomato", "lightsteelblue"), cex=1.5 )

dev.off()


png(paste(outdir,sampleID,"_segmentsgenomeplot.CNraw_loh.png",sep=""), width=1300,height=600)

par(mar = c(5, 5, 8, 4))

val=CNS$CN_raw
d=max(val)/100
plot(c(start,end),c(val,val),col="white",ylab="CN, CN Major, CN Minor",main=plottitle,ylim=c(-d,max(val)+4*d),xaxt='n',xlab="chromosomes")

for (i in 1:length(start)) {
  if (CNS$LOH[i]==1) {
    polygon(c(start[i],end[i],end[i],start[i]),c(min(val),min(val),max(val),max(val)),col="lightsteelblue",border="lightsteelblue",lwd=2)
  }
  
}

val=CNS$nBraw-d
segments(start,val,end,val,col="blue",lwd=4)
val=CNS$nAraw
segments(start,val,end,val,col="red",lwd=4)
val=CNS$CN_raw+d
segments(start,val,end,val,col="purple",lwd=4)
abline(v=y,col="grey")
posy=max(val)+2*d
i=1
l=0
for (x in chromo) {
  posx=(l+y[i])/2
  text(posx,posy,x,cex=1,srt=45)
  l=y[i]
  i=i+1
}
abline(h=CNS$ploidy[1],col="black",lty=2,lwd=1.5)

legend("topright", c("CN Total", "CN Major", "CN Minor", "LOH"),xpd=TRUE,horiz=T, inset=c(0,-0.1), bty = "n", lty=c(1,1,1,1), lwd=6, col = c("purple", "red", "blue", "lightsteelblue"), cex = 1)
dev.off()

