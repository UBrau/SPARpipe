#!/usr/bin/env Rscript
### Collect mappability stats files and output tables and graphs
###
### U. Braunschweig 2021

libMissing <- !require(optparse, quietly=T)
if (libMissing) {stop("Failed to load R package 'optparse'")}

#### Options
cArgs <- commandArgs(TRUE)

option.list <- list(
    make_option(c("-t", "--treatTab"),       action="store", type="character", metavar="FILE",
                help="Path of treatmtent file for experimental peaks. This file specifies all samples in the project,
                      including replicates. Required columns are ID, Replicate, Batch, Barcode.")
)
parser <- OptionParser(option_list=option.list,
                       usage="usage: %prog [options] OUTDIR",
                       description="Merge mappability from stats files, save tables and plots")
opt <- parse_args(parser, args=cArgs, positional_arguments=TRUE)

if (length(opt$args) == 0)
    stop(print_help(parser))
if(is.null(opt$options$treatTab)) {stop("Treatment table is required")}


### Check input
inDir <- sub("\\/*$","",opt$args[length(opt$args)]) # directories are checked upstream

if (!file.exists(opt$options$treatTab))            {stop("Treatment table not found")}
if (!dir.exists(inDir))                            {stop("Input directory not found")}
if (!dir.exists(file.path(inDir, "map")))          {stop("Subdirectory map/ not found")}
if (!dir.exists(file.path(inDir, "mappability")))  {stop("Subdirectory mappability/ not found")}

treat <- read.delim(opt$options$treatTab)


### Function defs

main <- function(workDir) {
    mapFilesF <- dir(file.path(workDir, "map"), pattern=".*(R1|Fwd|fwd|Fw|fw)\\.stats\\.txt")
    mapFilesR <- dir(file.path(workDir, "map"), pattern=".*(R2|Rev|rev|Rv|rv)\\.stats\\.txt")
    
    statsF <- lapply(file.path(workDir, "map", mapFilesF), read.delim, header=F)
    statsR <- lapply(file.path(workDir, "map", mapFilesR), read.delim, header=F)
    
    statsF <- do.call("rbind", statsF)
    statsR <- do.call("rbind", statsR)
    
    names(statsF) <- names(statsR) <- c("file","unmapped","total","mappability")
    rownames(statsF) <- rownames(statsR) <- sub("([^_]_W[0-9]+)_.+", "\\1", statsF[,1])
    
    totReads    <- rbind(statsF$total, statsR$total)
    mappability <- rbind(statsF$mappability, statsR$mappability)
    dimnames(totReads) <- dimnames(mappability) <- list(c("fwd","rev"), rownames(statsF))

    totReads    <- rbind(totReads, mean=colMeans(totReads))
    mappability <- rbind(mappability, mean=colMeans(mappability))
    mappedReads <- round(totReads * mappability / 100)

    write.csv(data.frame(Treatment=treat$Treatment, t(totReads)),
              file=file.path(workDir, "mappability", "TotalReads.csv"))
    write.csv(data.frame(Treatment=treat$Treatment, t(round(mappability, 2))),
              file=file.path(workDir, "mappability", "Mappability.csv"))
    write.csv(data.frame(Treatment=treat$Treatment, t(mappedReads)),
              file=file.path(workDir, "mappability", "MappedReads.csv"))

    
    ## Bar plots of mappability in order
    oriOrder <- breakLines(ncol(totReads), perLine=96)
            
    pdf(file.path(workDir, "mappability", "MappingStats.pdf"),
        wid=14, hei=4*length(oriOrder))
    par(mfrow=c(length(oriOrder),1), mar=c(8,4,3,3))

    ylim <- c(0, max(as.numeric(totReads[1:2,])))
    for (i in 1:length(oriOrder)) {
        barplot(totReads[1:2,oriOrder[[i]]], beside=T, las=3, col=c("grey40","grey80"), ylim=ylim,
                ylab="# reads", main="Raw reads",
                legend=T, args.legend=list(bty="n"))
    }
    
    ylim <- c(0, max(as.numeric(mappedReads[1:2,])))
    for (i in 1:length(oriOrder)) {
        barplot(mappedReads[1:2,oriOrder[[i]]], beside=T, las=3, col=c("grey40","grey80"), ylim=ylim,
                ylab="# reads", main="Mapped reads",
                legend=T, args.legend=list(bty="n"))
    }

    for (i in 1:length(oriOrder)) {    
        barplot(mappability[1:2,oriOrder[[i]]], beside=T, las=3, col=c("grey40","grey80"), ylim=c(0,100),
                ylab="% mappability", main="Mappability",
                legend=T, args.legend=list(bty="n"))
        abline(h=100, col="grey70")
    }
    dev.off()


    ## Scatter plot of reads vs mappability
    cols <- rep("black", nrow(treat))
    cols[treat$Type == "ctlNeg"] <- "dodgerblue"
    cols[treat$Type == "ctlPos"] <- "brown1"
    pch <- ifelse(grepl("ctl", treat$Type), 20, 1)
        
    pdf(file.path(workDir, "mappability", "MappabilityVsReads.pdf"),
        wid=6, hei=6.3)
    plot(totReads[3,], mappability[3,], ylim=c(0,100), type="n",
         xlab="Total Reads", ylab="% mappability", main="Mappability per sample")
    text(totReads[3,], mappability[3,], colnames(totReads), pos=1, cex=0.4, col="grey70")    
    points(totReads[3,], mappability[3,], cex=0.6, pch=pch, col=cols)    
    legend("topleft", pch=c(1,20,20,20), col=c("black","brown1","dodgerblue","black"),
           legend=c("Experimental","Pos. control","Neg. control","Other"), bty="n")
    dev.off()

    return(0)
}

breakLines <- function(x, perLine=96) {
    if (x <= perLine) {
        return(list(1:x))
    } else {
        nRows <- x %/% perLine + ceiling((x %% perLine) / x)
        se  <-  cbind(start = (0:(nRows - 1)) * perLine + 1,
                      end   = pmin(x, (1:nRows) * perLine)
                      )
       return(lapply(1:nrow(se), FUN=function(x) {se[x, "start"]:se[x, "end"]}))
    }
}



### Main

exit <- main(inDir)
if (exit != 0) {stop("Error when combining mappability stats")}
