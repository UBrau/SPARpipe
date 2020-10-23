#!/usr/bin/env Rscript
### Combine read counts and RPM from raw counts per well for SARS-CoV-2 project
### U. Braunschweig 2020


libMissing <- !require(optparse, quietly=T)
if (libMissing) {stop("Failed to load R package 'optparse'")}

#### Options
args <- commandArgs(TRUE)

option.list <- list(
    make_option(c("-e", "--eventTab"), type="character",
        help="CSV file describing the amplicons - see README"),
    make_option(c("-t", "--treatTab"), type="character",
        help="CSV file describing the samples - see README"),
    make_option(c("-d", "--direction"),  default="both",
                help="Which reads to use. 'fwd' for forward/R1, 'rev' for reverse/R1, or 'both for mean [%default]"),
    make_option(c("-i", "--imbal"),  default=1.5,
                help="Maximum imbalance between fwd and rev reads before a warning is printed [%default]"),
    make_option(c("-b", "--imbalCts"),  default=20,
                help="Minimum reads (average of fwd and rev) to do imbalance test [%default]"),
    make_option(c("-c", "--cores"),  default=1,
                help="Number of CPUs [%default]")

)
parser <- OptionParser(option_list=option.list,
                       usage="usage: %prog [options] OUTDIR",
                       description="Combine raw reads from invidual samples for expression analysis")
opt <- parse_args(parser, args=args, positional_arguments=1)

inDir <- sub("/*$","", opt$args[1])

libMissing <- !require(parallel, quietly=T)
if (libMissing) {stop("Failed to load R package 'parallel'")}


### Function defs

combineCounts <- function(inDir, treat, junc, cores=1) {
### Reads individual raw well counts from the 'counts' folder and combine into one table
    
    files <- list.files(file.path(inDir, "counts"), pattern=paste(treat$Batch[1], ".*W.*.counts.tab", sep=""), full.names=T)
    fileBC <- sub(".*(W[0-9]+).*", "\\1", files)
    
    if (length(files) == 0)                {stop("Raw well data not found")}
    if (!(all(treat$Barcode %in% fileBC))) {stop("Not all raw well data files were found")}
    if (!all(fileBC == treat$Barcode))     {stop("Mismatch betwetween treatment table and found raw well data files")}
    
    ## Get total reads
    mapFilesF <- dir(file.path(inDir, "map"), pattern=".*(R1|Fwd|fwd|Fw|fw)\\.stats\\.txt")
    mapFilesR <- dir(file.path(inDir, "map"), pattern=".*(Rw|Rev|rev|Rv|rv)\\.stats\\.txt")    
    statsF <- lapply(file.path(inDir, "map", mapFilesF), read.delim, header=F)
    statsR <- lapply(file.path(inDir, "map", mapFilesR), read.delim, header=F)
    statsF <- do.call("rbind", statsF)
    statsR <- do.call("rbind", statsR)

    bcF <- sub(".+(W[0-9]+).*", "\\1", statsF[,1])
    bcR <- sub(".+(W[0-9]+).*", "\\1", statsR[,1])
    if (!all(bcF == bcR))             {stop("Mismatch between fwd and rev read stats files")}
    if (!all(treat$Barcode %in% bcF)) {stop("Not all expected mapping stats files were found")}    
    names(statsF) <- names(statsR) <- c("file","unmapped","total","mappability")
    rownames(statsF) <- rownames(statsR) <- sub("[^_]+_(.+)_[^_]+", "\\1", statsF[,1])
    statsF <- statsF[match(treat$Barcode, bcF),]
    statsR <- statsF[match(treat$Barcode, bcR),]
    
    input  <- mclapply(files, read.csv, sep="\t", header=F, mc.cores=cores)
    counts <- sapply(input, "[[", 4)
    inputLabel <- paste(input[[1]][,1], input[[1]][,2], sep="_")
    geneCounts <- aggregate(counts, by=list(Gene=paste(input[[1]][,1], substr(input[[1]][,2], 1, 2), sep=".")), FUN=sum)
    
    if (any(!is.na(junc$Label))) {
        collapseJunc <- junc[!is.na(junc$Label),]
        juncF <- strsplit(collapseJunc$JunctionsFw, split=",")
        juncR <- strsplit(collapseJunc$JunctionsRv, split=",")
        
        collapseLabel <- rep(NA, nrow(input[[1]]))
        for (i in 1:nrow(collapseJunc)) {
            collapseLabel[input[[1]][,1] == collapseJunc$Gene[i] & input[[1]][,2] %in% juncF[[i]]] <-
                paste0(collapseJunc$Gene[i], ".fw.", collapseJunc$Label[i])
            collapseLabel[input[[1]][,1] == collapseJunc$Gene[i] & input[[1]][,2] %in% juncR[[i]]] <-
                paste0(collapseJunc$Gene[i], ".rv.", collapseJunc$Label[i])
        }
        eventCounts <- aggregate(counts, by=list(Event=collapseLabel), FUN=sum)
        eventCounts <- data.frame(Gene=sub("(.*)\\.[^.]+", "\\1", eventCounts$Event), eventCounts)
        subtractCounts <- aggregate(eventCounts[,-c(1,2)], by=list(Gene=eventCounts$Gene), FUN=sum)
        
        out <- geneCounts[!(geneCounts$Gene %in% subtractCounts$Gene),]
        tmp <- geneCounts[match(subtractCounts$Gene, geneCounts$Gene),]
        tmp <- data.frame(Gene=paste0(subtractCounts$Gene, ".RNA"),
                          as.matrix(tmp[,-1]) - as.matrix(subtractCounts[,-1])
                          )
        out <- rbind(out,
                     tmp,
                     data.frame(Gene=sub("A1","GEN",eventCounts$Event), eventCounts[,-c(1,2)])
                     )
        out <- out[order(out$Gene),]
    } else {
        out <- geneCounts
    }
    outNames <- as.character(out$Gene)
    out <- t(out[,-1])
    colnames(out) <- outNames
    rownames(out) <- treat$Barcode

    out <- cbind(out,
                 totReads.fw    = statsF[,"total"],
                 totReads.rv    = statsR[,"total"],
                 mappedReads.fw = statsF[,"total"] - statsF["unmapped"],
                 mappedReads.rv = statsR[,"total"] - statsR["unmapped"],
                 mappability.fw = round(100 * (1 - (statsF[,"unmapped"] / statsF[,"total"])), 2),
                 mappability.rv = round(100 * (1 - (statsR[,"unmapped"] / statsR[,"total"])), 2)
                 )

    out
}

avgFRcounts <- function(x, junc, direction=c("both","fwd","rev")[1], imbal=1.2, imbalCts=10) {
### Average forward and reverse read counts, or select one; Also check for imbalance more than 'imbal'
    eventsFw <- sort(grep(".*\\.fw", colnames(x), value=T))
    fwInd <- sapply(eventsFw, FUN=function(y) {grep(y, colnames(x))})
    rvInd <- sapply(sub("\\.fw", ".rv", eventsFw), FUN=function(y) {grep(y, colnames(x))})

    ## Check balance of fwd and rev
    bal <- abs(log(10 + x[,fwInd]) / log(10 + x[,rvInd]))
    meanFR <- (x[,fwInd] + x[,rvInd]) / 2
    warn <- (bal < 1 - log(imbal) | bal > 1 + log(imbal)) & meanFR > imbalCts
    if (any(as.vector(warn))) {
        for (i in which(apply(warn, MAR=2, any))) {
            warning(length(which(warn[,i])), " imbalanced samples for ", colnames(meanFR)[i])
        }
    }    

    if (direction == "fwd") {
        out <- x[,fwInd]
        colnames(out) <- sub("\\.fw", "", colnames(out))
    }
    if (direction == "rev") {
        out <- x[,rvInd]
        colnames(out) <- sub("\\.rv", "", colnames(out))
    }
    if (direction == "both") {
        if (is.list(fwInd) || is.list(rvInd)) {
            stop("Forward and reverse read columns in raw counts file could not be matched for all events")
        }
        
        out <- round(meanFR)
        colnames(out) <- sub("\\.fw", "", colnames(out))
    }

    out
}

main <- function(inDir, treat, junc, direction=c("both","fwd","rev")[1], cores=1, imbal=1.2, imbalCts=10) {
### Collect mapped read numbers from individual samples and creat tables with
### raw fwd and rev counts, average counts (if selected), and RPM

    rawCounts <- combineCounts(inDir, treat, junc, cores=cores)
    counts <- avgFRcounts(rawCounts, direction=direction, imbal=imbal,
                          junc = junc)
    tmp <- counts[,-which(colnames(counts) %in% c("mappability","totReads"))]
    rpm <- round(1000000 * tmp / rowSums(tmp), 1)

    treatCols <- unlist(sapply(c("ID","Barcode","Batch","Replicate","Treatment","Type","Status"),
                               FUN=grep, x=names(treat)))
    treat <- treat[,treatCols]

    if (is.null(treat$Batch)) {
        outName <- ""
    } else {
        outName <- paste0(treat$Batch[1], "_") 
    }
    dirTag <- paste0("_", direction)

    outDir <- file.path(inDir, "expression")
    if (!dir.exists(outDir)) {dir.create(outDir, recursive=TRUE)}
    
    write.csv(data.frame(treat, rawCounts, check.names=F),
              file=file.path(outDir, paste0(outName, "CountsRaw.csv")))
    write.csv(data.frame(treat, counts, check.names=F),
              file=file.path(outDir, paste0(outName, "Counts", dirTag, ".csv")), row.names=F)
    write.csv(data.frame(treat, rpm, check.names=F),
              file=file.path(outDir, paste0(outName, "RPM", dirTag, ".csv")), row.names=F)

    cat("Read and RPM files saved in", outDir, "\n")
}





### Body

## Check input
if (length(opt$args) != 1)
    stop(print_help(parser))
if (is.null(opt$options$eventTab)) 
    stop("eventTab must be provided")
if (is.null(opt$options$treatTab)) 
    stop("treatTab must be provided")
if (!(opt$options$direction %in% c("both","fwd","rev"))) 
    stop("direction must be one of 'both', 'fwd', or 'rev'")
    
## Load tables
treat <- read.delim(opt$options$treatTab)
treat <- treat[order(treat$Barcode),]

junc <- read.delim(opt$options$eventTab, as.is=T)


## Body
main(inDir, treat, junc,
     direction = opt$options$direction,
     cores     = opt$options$cores,
     imbal     = opt$options$imbal,
     imbalCts  = opt$options$imbal
     )

#opt <- list(
#    options = list(
#        treatTab  = "~/bin/SPARpipe/examples/COVID-19/SampleTable.tab",
#        eventTab  = "~/bin/SPARpipe/examples/COVID-19/amplicons/EventJunctions.tab",
#        direction = "fwd",
#        imbal     = 1.5,
#        imbalCts  = 20,
#        cores     = 1
#    ),
#    args    = "~/bin/SPARpipe/examples/COVID-19/test_align/"
#)

