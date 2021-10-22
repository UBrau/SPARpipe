#!/usr/bin/env Rscript
### Combine PSI and RPM values from all experiments into one super file
###
### K. Ha & U. Braunschweig 2013-2020
### Changes: - More intuitively named files and columns
###          - Adapt to different file structure of input
###          - Fix issue with simultaneous processing of batches


libMissing <- !require(optparse, quietly=T)
libMissing <- libMissing | !require(parallel, quietly=T)
libMissing <- libMissing | !require(stringr, quietly=T)
libMissing <- libMissing | !require(plyr, quietly=T)
if (libMissing) {stop("Failed to load R package(s) in merge_tables.R")}

### Parse input
args <- commandArgs(TRUE)
option.list <- list(
    make_option(c("-b", "--batch"),  default="",
        help="Batch identifier [%default]"),
    make_option(c("-c", "--cores"),  default=1,
        help="Number of CPUs [%default]")
)
parser <- OptionParser(option_list=option.list,
                       usage="usage: %prog [options] OUTDIR",
                       description="Merge PSI and RPM table from multiple samples from the same batch")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) == 0)
    stop(print_help(parser))

inDir <- sub("/*$","",opt$args[1])


### Function definitions
main <- function(batch, cores=1) {
    psi.files <- list.files(file.path(inDir, "welldata"), pattern=paste(batch, ".*W.*.RAW.tab", sep=""), full.names=T)

    psi.input     <- mclapply(psi.files, read.csv, sep="\t", mc.cores=cores)
    rpm.input     <- mclapply(psi.input, "[", c("Event","Element", "RPM.Fw", "RPM.Rv", "RPM"),
                              mc.cores=cores)
    read.input    <- mclapply(psi.input, "[", c(which(names(psi.input[[1]]) == "Event"),
                                                which(names(psi.input[[1]]) == "Element"),
                                                grep("Reads", names(psi.input[[1]]))), 
                              mc.cores=cores)
    countIE.input <- mclapply(psi.input, "[", c("Event", "Element", "Counts.InEx"),
                              mc.cores=cores)
    psi.input     <- mclapply(psi.input, "[", c("Event", "Element", "PSI.Fw", "SD.Fw", "PSI.Rv", "SD.Rv", "PSI"),
                              mc.cores=cores)
    names(psi.input) <- names(read.input) <- names(countIE.input) <- names(rpm.input) <- sub(".*_(W[0-9]+).+", "\\1", psi.files)

    targetID <- paste(psi.input[[1]]$Event, psi.input[[1]]$Element, sep=".")


    ## PSI
    super.psi <- matrix(nrow = nrow(psi.input[[1]]), ncol=5 * length(psi.input),
                        dimnames=list(c(), paste(rep(names(psi.input), each=5), names(psi.input[[1]])[3:7], sep="."))
                        )
    for (i in 1:length(psi.input)) {
        super.psi[,(i - 1) * 5 + 1:5] <- as.matrix(psi.input[[i]][,3:7])
    }
    super.psi <- data.frame(Event = targetID, super.psi)

    write.table(super.psi, file=file.path(inDir, paste("batchdata/PSI.FULL_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    write.table(super.psi[,c(1,seq(6,ncol(super.psi),5))], 
        file=file.path(inDir, paste("batchdata/PSI_", batch, ".tab", sep="")),
        sep="\t", quote=F, row.names=F)
    cat("Done merging PSI files\n")

    
    ## RPM
    super.rpm <- matrix(nrow = nrow(psi.input[[1]]), ncol=3 * length(psi.input),
                        dimnames=list(c(), paste(rep(names(rpm.input), each=3), names(rpm.input[[1]])[3:5], sep="."))
                        )
    for (i in 1:length(rpm.input)) {
        super.rpm[,(i - 1) * 3 + 1:3] <- as.matrix(rpm.input[[i]][,3:5])
    }
    super.rpm <- data.frame(Event = targetID, super.rpm)

    write.table(super.rpm, file=file.path(inDir, paste("batchdata/RPM.FULL_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    write.table(super.rpm[,c(1,seq(4,ncol(super.rpm),3))], 
                file=file.path(inDir, paste("batchdata/RPM_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    cat("Done merging RPM files\n")


    ## Reads
    super.reads <- matrix(nrow = nrow(psi.input[[1]]), ncol=2 * length(psi.input),
                        dimnames=list(c(), paste(rep(names(read.input), each=2), names(read.input[[1]])[3:4], sep="."))
                        )
    for (i in 1:length(read.input)) {
        super.reads[,(i - 1) * 2 + 1:2] <- as.matrix(read.input[[i]][,3:4])
    }
    super.reads <- data.frame(Event = targetID, super.reads)

    write.table(super.reads, file=file.path(inDir, paste("batchdata/ReadsPerEvent_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    cat("Done merging reads files\n")


    ## Pseudo-inclusion/exclusion reads
    super.ie <- matrix(nrow = nrow(psi.input[[1]]), ncol=length(psi.input),
                       dimnames=list(c(),
                                     paste(names(read.input), names(countIE.input[[1]])[3], sep=".")
                                     )
                        )
    for (i in 1:length(read.input)) {
        super.ie[,i] <- as.matrix(countIE.input[[i]][,3])
    }
    super.ie <- data.frame(Event = targetID, super.reads)

    write.table(super.ie, file=file.path(inDir, paste("batchdata/InclExclCounts_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    cat("Done merging inclusion/exclusion read count files\n")

    return(0)
}

### Meat
exit <- main(opt$options$batch, cores=opt$options$cores)
if (exit != 0) {stop("Error during merging")}
