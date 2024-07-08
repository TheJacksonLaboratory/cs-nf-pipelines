#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#1Mb windows, sliding every 0.5Mb
file <- get(load(args[1]))

#get windows with NA

win=data.frame(matrix(ncol=0,nrow=0))

for (i in 1:length(file$ratio)) {
    tmp=file$ratio[[i]][file$ratio[[i]]$N==0,]
    if (nrow(tmp)>0) {
        if (i==23) i="X"
        else if (i==24) i="Y"
        win1=cbind(paste("chr",i,sep=""),tmp)
        colnames(win1)[1] <- "chromo"
        win <- rbind(win,win1)
    }
}

write.table(win,file=args[2], sep="\t",row.names=F,col.names=F,quote=F)
