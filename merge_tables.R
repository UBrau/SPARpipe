#!/usr/bin/env Rscript
### Combine PSI and RPM values from all experiments into one super file
###
### K. Ha 2013-11-14
### Modified U. Braunschweig 10/2014 and 04/2017
### Changes: - More intuitively named files and columns
###          - Adapt to different file structure of input


libMissing <- !require(optparse, quietly=T)
libMissing <- libMissing | !require(parallel, quietly=T)
libMissing <- libMissing | !require(stringr, quietly=T)
libMissing <- libMissing | !require(plyr, quietly=T)
if (libMissing) {stop("Failed to load R package(s) in merge_tables.R")}

### Parse input
args <- commandArgs(TRUE)
option.list <- list(
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
main <- function(batch) {
    psi.files <- list.files(file.path(inDir, "welldata"), pattern="W.*.RAW.tab", full.names=T)

    psi.input     <- lapply(psi.files, read.csv, sep="\t")
    rpm.input     <- lapply(psi.input, "[", c("Event", "RPM.Fw", "RPM.Rv", "RPM"))
    read.input    <- lapply(psi.input, "[", c(which(names(psi.input[[1]]) == "Event"), grep("Reads", names(psi.input[[1]]))))
    countIE.input <- lapply(psi.input, "[", c("Event", "Counts.InEx"))
    psi.input     <- lapply(psi.input, "[", c("Event", "PSI.Fw", "SD.Fw", "PSI.Rv", "SD.Rv", "PSI"))
    names(psi.input) <- names(read.input) <- names(countIE.input) <- names(rpm.input) <- str_extract(psi.files, "W\\d+")

    ## PSI
    super.psi <- join_all(psi.input, by="Event")
    colnames(super.psi)[2:ncol(super.psi)] <- paste(rep(names(psi.input), each=ncol(psi.input[[1]])-1),
                                                    colnames(super.psi)[2:ncol(super.psi)],
                                                    sep=".")
    write.table(super.psi, file=file.path(inDir, paste("batchdata/PSI.FULL_", batch, ".tab", sep="")), sep="\t", quote=F, row.names=F)
    write.table(super.psi[,c(1,seq(6,ncol(super.psi),5))], 
        file=file.path(inDir, paste("batchdata/PSI_", batch, ".tab", sep="")),
        sep="\t", quote=F, row.names=F)
    cat("Done merging PSI files\n")


    ## RPM
    super.rpm <- join_all(rpm.input, by="Event")

    colnames(super.rpm)[2:ncol(super.rpm)] <- paste(rep(names(rpm.input), 
                                                        each=ncol(rpm.input[[1]])-1),
                                                    colnames(super.rpm)[2:ncol(super.rpm)],
                                                    sep=".")

    super.rpm <- data.frame(Event=super.rpm$Event, super.rpm[,2:ncol(super.rpm)])

    write.table(super.rpm, file=file.path(inDir, paste("batchdata/RPM.FULL_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    write.table(super.rpm[,c(1,seq(5,ncol(super.rpm),3))], 
                file=file.path(inDir, paste("batchdata/RPM_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    cat("Done merging RPM files\n")


    ## Reads
    super.reads <- join_all(read.input, by="Event")

    colnames(super.reads)[2:ncol(super.reads)] <- paste(rep(names(read.input), 
                                                        each=ncol(read.input[[1]])-1),
                                                    colnames(super.reads)[2:ncol(super.reads)],
                                                    sep=".")

    write.table(super.reads, file=file.path(inDir, paste("batchdata/ReadsPerEvent_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    cat("Done merging reads files\n")


    ## Pseudo-inclusion/exclusion reads
    super.ie <- join_all(countIE.input, by="Event")

    colnames(super.ie)[2:ncol(super.ie)] <- paste(rep(names(countIE.input), 
                                                        each=ncol(countIE.input[[1]])-1),
                                                    colnames(super.ie)[2:ncol(super.ie)],
                                                    sep=".")

    write.table(super.ie, file=file.path(inDir, paste("batchdata/InclExclCounts_", batch, ".tab", sep="")),
                sep="\t", quote=F, row.names=F)
    cat("Done merging inclusion/exclusion read count files\n")

    return(0)
}

### Meat
batch <- sort(unique(sub("_$", "", sub("(.*)W[0-9].+", "\\1", dir(file.path(inDir, "welldata"), pattern=".*W[0-9]+.+RAW.tab")))))
exit <- mclapply(batch, main, mc.cores=min(opt$options$cores, length(batch)))
if (!all(exit == 0)) {stop("Error during merging")}