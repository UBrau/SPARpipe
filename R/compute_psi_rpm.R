#!/usr/bin/env Rscript
### Tabulate coverage of multiplex samples and compute PSI and RPM
###
### K. Ha 2013-11-08
### Modified U. Braunschweig 10/2014 and 04/2017
### Changes: - Removed calculation of legacy quality values
###          - 'up' and 'dn' renamed to 'fw' and 'rv' throughout
###          - changed order of columns in output
###          - allow parallel processing


libMissing <- !require(optparse, quietly=T)
if (libMissing) {stop("Failed to load R package 'optparse'")}

#### Options
args <- commandArgs(FALSE)

option.list <- list(
    make_option(c("-e", "--events"), type="character",
        default="/home/blencowe/blencowe31/ulrich/proj/SPARseq.HongAndy.neural.170202/seq/junctions/Events108_11-101_UB170417.tab",
        help="AS event junction table [%default]"),
    make_option(c("-c", "--cores"),  default=1,
        help="Number of CPUs [%default]")
)
parser <- OptionParser(option_list=option.list,
                       usage="usage: %prog [options] OUTDIR",
                       description="Tabulate coverage of multiplex samples to compute PSI and RPM")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) == 0)
    stop(print_help(parser))

libMissing <- !require(parallel, quietly=T)
if (libMissing) {stop("Failed to load R package 'parallel'")}


### Get gene AS type and junctions to use
AS.types <- read.table(opt$options$events, sep="\t", header=T)

### Functions
get_beta_sd <- function(incl, excl) {
### Calculate SD from alpha and beta of beta distr.
    incl <- incl + 1
    excl <- excl + 1
    s <- (incl*excl) / ((incl+excl)^2*(incl+excl+1))
    return(sqrt(s))
}

get_counts.IE <- function(readsUp, readsDn, psi) {
### Get pseudo shape parameters of a beta distribution based on PSI and the number of reads
### of either the greater of the fwd and rev reads or, if only one is used, that one.
    useUp <- readsUp >= readsDn
    useUp[is.na(useUp) & !is.na(readsUp)] <- TRUE
    useUp[is.na(useUp) & !is.na(readsDn)] <- FALSE
    reads <- ifelse(useUp, readsUp, readsDn)

    alpha <- reads * psi
    beta  <- reads * (1 - psi)

    paste(ifelse(is.na(alpha), NA, round(alpha, 2)), 
          ifelse(is.na(beta),  NA, round(beta, 2)),
	  sep="=")      
}

compute_rpm <- function(m) {
    cUp <- grep("Counts\\.fw", names(m), value=T)
    cDn <- grep("Counts\\.rv", names(m), value=T)
    total.up <- sum(m[,cUp], na.rm=T)
    total.dn <- sum(m[,cDn], na.rm=T)
    
    newm <- m
    newm[,c("Reads.Fw", "Reads.Rv", "RPM.Fw", "RPM.Rv", "RPM")] <- NA
    
    newm$Reads.Fw <- rowSums(newm[,cUp], na.rm=T)
    newm$Reads.Rv <- rowSums(newm[,cDn], na.rm=T)
    newm <- within(newm,
                   {
                       RPM.Fw <- 1000000 * Reads.Fw / total.up
                       RPM.Rv <- 1000000 * Reads.Rv / total.dn
                       RPM <- 1000000 * (Reads.Fw + Reads.Rv) / (total.up + total.dn)
                   })
    
    newm <- within(newm,
                   {
                       RPM.Fw <- round(RPM.Fw, 2)
                       RPM.Rv <- round(RPM.Rv, 2)
                       RPM <- round(RPM, 2)
                   })
    return(newm)
}

compute_psi <- function(m) {  
    newm <- m
    newm[,c("PSI.Fw", "SD.Fw","RPM.Fw","Reads.Fw","PSI.Rv", "SD.Rv", "RPM.Rv","Reads.Rv","PSI","Counts.InEx","RPM")] <- NA

    upJunctions <- strsplit(as.character(newm$JunctionsFw), split=",")
    dnJunctions <- strsplit(as.character(newm$JunctionsRv), split=",")
    upJunctions <- ifelse(sapply(upJunctions, FUN=function(x) {is.na(x[1])}),
                          yes = NA,
                          no  = lapply(upJunctions, FUN=function(x) {paste("Counts", x, sep=".")})
                          )
    dnJunctions <- ifelse(sapply(dnJunctions, FUN=function(x) {is.na(x[1])}),
                          yes = NA,
                          no  = lapply(dnJunctions, FUN=function(x) {paste("Counts", x, sep=".")})
                          )
    allUpJunctions <- grep("Counts\\.fw", names(newm), value=T)
    allDnJunctions <- grep("Counts\\.rv", names(newm), value=T)

    ## count inclusion and exclusion reads
    inclUp <- sapply(1:nrow(newm), FUN=function(x) {
        if (is.na(upJunctions[[x]][1])) {
            return(NA)
        } else {
            return(sum(newm[x,upJunctions[[x]]]))
        }
    })
    inclDn <- sapply(1:nrow(newm), FUN=function(x) {
        if (is.na(dnJunctions[[x]][1])) {
            return(NA)
        } else {
            return(sum(newm[x,dnJunctions[[x]]]))
        }
    })
    exclUp <- sapply(1:nrow(newm), FUN=function(x) {
        if (is.na(upJunctions[[x]][1])) {
            return(NA)
        } else {
            return(sum(newm[x,setdiff(allUpJunctions, upJunctions[[x]])], na.rm=T))
        }
    })
    exclDn <- sapply(1:nrow(newm), FUN=function(x) {
        if (is.na(dnJunctions[[x]][1])) {
            return(NA)
        } else {
            return(sum(newm[x,setdiff(allDnJunctions, dnJunctions[[x]])], na.rm=T))
        }
    })
    
    newm$PSI.Fw      <- inclUp / (inclUp + exclUp)
    newm$PSI.Rv      <- inclDn / (inclDn + exclDn)
    newm$PSI         <- rowMeans(newm[,c("PSI.Fw","PSI.Rv")], na.rm=T)
    newm$SD.Fw       <- get_beta_sd(inclUp, exclUp)
    newm$SD.Rv       <- get_beta_sd(inclDn, exclDn)
    newm$Counts.InEx <- get_counts.IE(inclUp + exclUp, inclDn + exclDn, newm$PSI)
    
    newm <- within(newm,
                   {
                       PSI.Fw <- round(PSI.Fw * 100, 2)
                       PSI.Rv <- round(PSI.Rv * 100, 2)
                       PSI    <- round(PSI    * 100, 2)
                       PSI[is.na(PSI)] <- NA
                   })
    
    return(newm)
}

main <- function(file) {
    m <- read.table(file, sep="\t", col.names=c("Gene", "Junction", "Label", "Counts"))
    m <- m[,c(1,2,4)]
    newm <- reshape(m, v.names="Counts", idvar="Gene", timevar="Junction",
                    direction="wide")
    newm <- merge(AS.types, newm, by=1)  # duplicate lines with multiple events
    newm$Event <- ifelse(is.na(newm$Event), yes=as.character(newm$Gene), no=paste(newm$Gene, newm$Event, sep="."))
    newm <- compute_psi(newm)
    newm <- compute_rpm(newm)
    
    newm <- within(newm,
                   {
                       PSI.Fw <- sprintf("%.2f", PSI.Fw)
                       PSI.Rv <- sprintf("%.2f", PSI.Rv)
                       PSI    <- sprintf("%.2f", PSI)
                       RPM.Fw <- sprintf("%.2f", RPM.Fw)
                       RPM.Rv <- sprintf("%.2f", RPM.Rv)
                       SD.Fw  <- sprintf("%.2f", SD.Fw)
                       SD.Rv  <- sprintf("%.2f", SD.Rv)
                       RPM    <- sprintf("%.2f", RPM)
                   })
    outfile <- file.path(outDir, "welldata", paste(sub("(.*W[0-9]+[^.]+).+", "\\1", basename(file)), "RAW.tab", sep="."))
    write.table(newm, file=outfile, sep="\t", quote=F, row.names=F)
    cat("Converted", file, "to", outfile, "\n")
}


### Check input
outDir  <- sub("/*$", "", opt$args[length(opt$args)])

if (!dir.exists(file.path(outDir, "counts"))) {
    stop("Directory ", outDir, "/counts not found")
}                                                                                                                                                                           

files <- dir(file.path(outDir, "counts"), pattern=".*W[0-9]+.+counts.tab")
if (length(files) < 1) {
   stop("No counts files found in ", outDir, "/counts")
}  


### Execute commands
r <- mclapply(file.path(outDir, "counts", files), main, mc.cores=opt$options$cores)
