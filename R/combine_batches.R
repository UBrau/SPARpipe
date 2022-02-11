#!/usr/bin/env Rscript
### Part of the SPAR-seq pipeline
### After generating raw PSI reads and shape counts, combine the different
### 'batches' (e.g., lanes), re-calculate the PSI of Alt3/Alt5 as a fraction of
### the host exon if applicable, perform a normalization, and calculate
### average PSI, dPSI, and SSMD for each treatment.

cArgs <- commandArgs(TRUE)

## Check dependencies
libMissing <- !require("optparse", quietly=T) && stop("Failed to load R package 'optparse")


### Capture input
opt.list <- list(
    make_option(c("-e", "--eventTab"),       action="store", type="character", metavar="FILE",
                help="Names of event table, tab-separated."),
    make_option(c("-t", "--treatTab"),       action="store", type="character", metavar="FILE",
                help="Path of treatmtent file for experimental peaks. This file specifies all samples in the project,
                      including replicates. Required columns are ID, Replicate, Batch, Barcode.")
)
parser <- OptionParser(option_list=opt.list,
                       usage="usage: %prog [options] OUTDIR",
                       description="Merge PSI, RPM, inclusion-exclusion read, and reads per event table from all batches")
opt <- parse_args(parser, args=cArgs, positional_arguments=TRUE)

if (length(opt$args) == 0)
    stop(print_help(parser))
if(is.null(opt$options$eventTab)) {stop("Event table is required")}
if(is.null(opt$options$treatTab)) {stop("Treatment table is required")}


### Function definitions
getTable <- function(x, ev, files) {
### Read in PSI tables per batch and reformat
### x: a line number reffering to the files table
    raw <- read.delim(file.path(inDir, "batchdata", files$name[x]), as.is=T)
    if (!all(raw$Event == ev$Name)) {
        stop("Events in ", files$name[x], " not identical with those in events table")
    }
    out <- data.frame(Batch   = files$batch[x],
                      Barcode = sub("(W[0-9]+).*", "\\1", names(raw)[-1]),
                      t(raw[,-1]),
                      stringsAsFactors=F
                      )
    if (files$type[x] == "RPE") {out <- out[seq(1, nrow(out) - 1, 2),]}
    names(out)[-c(1,2)] <- as.character(raw$Event)
    out
}

cTables <- function(x, type="", treat) {
### Concatenate a list of tables
### x: a List of reformatted tables
    if (length(x) == 1) {
        out <- x[[1]]
    } else {
        out <- x[[1]]
        for (i in 2:length(x)) {out <- rbind(out, x[[i]])}
    }
    if (nrow(treat) != nrow(out) || !all(paste(treat$Batch, treat$Barcode) %in% paste(out$Batch, out$Barcode))) {
        stop("Rows in combined ", type, " table (n=", nrow(out), 
            ") do not match rows in treatments table (n=", nrow(treat), "). Check treatTab!")
    }
    key <- merge(data.frame(Sample = paste(treat$Batch, treat$Barcode, sep="."), treatInd = 1:nrow(treat)),
                 data.frame(Sample = paste(out$Batch, out$Barcode, sep="."), outInd = 1:nrow(out)),
                 by.x=1, by.y=1)
    key <- key[order(key$treatInd),]
    out <- out[key$outInd,]
    data.frame(Replicate = treat$Replicate,
               ID        = treat$ID,
               out[,-c(1,2)],
               check.names=F
               )
}

depEvents <- function(x, ev) {
### Re-normalize dependent events by the PSI of their parent
    rel <- which(!is.na(ev$Relative))
    dep <- sapply(ev$Name[rel], grep, x=names(x))
    if (length(dep) > 0)  {
        parName <- ifelse(is.na(ev$Relative), as.character(ev$Event), paste(ev$Event, ev$Relative, sep="."))
        if (!all(parName[rel] %in% names(x))) {
            warning("Parent event not found for ", 
	        paste(parName[!is.na(ev$Relative) & !(parName %in% names(x))], collapse=", ")
	    )
	    dep <- dep[parName[rel] %in% names(x)]
	    rel <- rel[parName[rel] %in% names(x)]
        }
        par <- sapply(parName[rel], grep, x=names(x))
        for (i in 1:length(dep)) {
            new <- round(x[,dep[i]] * x[,par[i]] / 100, 2)
            if (any((new < 0 | new > 100) & !is.na(new)) & !all(is.na(new))) {
                warning("Relative PSI outside range for event ", names(x)[dep[i]], ": ", min(new, na.rm=T),
                        " - ", max(new, na.rm=T))
            }
            new[new < 0]    <-   0
            new[new > 100]  <- 100
            x[,dep[i]] <- new
        }    
    }
    return(x)
}


### Check input
inDir <- sub("\\/*$","",opt$args[length(opt$args)]) # directories are checked upstream

if (!file.exists(opt$options$eventTab))         {stop("Events table not found")}
if (!file.exists(opt$options$treatTab))         {stop("Treatments table not found")}
if (!dir.exists(inDir))                         {stop("Input directory not found")}
if (!dir.exists(file.path(inDir, "batchdata"))) {stop("Subdirectory batchdata/ not found")}

batches <- sort(unique(sub("PSI_(.+)\\.tab", "\\1", dir(file.path(inDir, "batchdata"), pattern="PSI_.*tab"))))
cat(length(batches), "batch(es) found:", paste(batches, collapse=", "), "\n")

psiFiles       <- data.frame(name  = paste("PSI_", batches, ".tab", sep=""),
                             type  = "PSI",
                             batch = batches)
psiFiles$found <- psiFiles$name %in% dir(file.path(inDir, "batchdata"))
rpeFiles       <- data.frame(name  = paste("ReadsPerEvent_", batches, ".tab", sep=""),
                             type  = "RPE",
                             batch = batches)
rpeFiles$found <- rpeFiles$name %in% dir(file.path(inDir, "batchdata"))
iecFiles       <- data.frame(name  = paste("InclExclCounts_", batches, ".tab", sep=""),
                             type  = "IEC",
                             batch = batches)
iecFiles$found <- iecFiles$name %in% dir(file.path(inDir, "batchdata"))
rpmFiles       <- data.frame(name  = paste("RPM_", batches, ".tab", sep=""),
                             type  = "RPM",
                             batch = batches)
rpmFiles$found <- rpmFiles$name %in% dir(file.path(inDir, "batchdata"))
files <- rbind(psiFiles, rpeFiles, iecFiles, rpmFiles)

if (any(!files$found)) {
    stop("Batch data file(s) not found: ", paste(files$name[!files$found], collapse=", "))
}


### Load, combine and save tables
treat <- read.delim(opt$options$treatTab)
ev    <- read.delim(opt$options$eventTab, sep="\t", header=T, comment.char="#")


ev$Name <- ifelse(is.na(ev$Element), as.character(ev$Event), paste(ev$Event, ev$Element, sep="."))
ev <- ev[order(ev$Name),]

psi <- lapply(which(files$type == "PSI"), getTable, ev=ev, files=files)
psi <- cTables(psi, type="PSI", treat=treat)
psi <- depEvents(psi, ev=ev)
write.table(psi, row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "raw", "PSI.raw.tab"))

rpe <- lapply(which(files$type == "RPE"), getTable, ev=ev, files=files)
rpe <- cTables(rpe, type="RPE", treat=treat)
write.table(rpe, row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "raw", "ReadsPerEvent.tab"))

iec <- lapply(which(files$type == "IEC"), getTable, ev=ev, files=files)
iec <- cTables(iec, type="IEC", treat=treat)
write.table(iec, row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "raw", "InclExclCounts.tab"))

rpm <- lapply(which(files$type == "RPM"), getTable, ev=ev, files=files)
rpm <- cTables(rpm, type="RPM", treat=treat)
write.table(rpm, row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "raw", "RPM.tab"))

