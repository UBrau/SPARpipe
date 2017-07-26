#!/usr/bin/env Rscript
### Normalize raw PSI and compute SSMD
### Part of the SPARpipe pipeline for analysis of SPAR-seq data
###
### U. Braunschweig, 07/2017


libMissing <- !require(optparse, quietly=T) && stop("Failed to load R package 'optparse'")

### Parse input
args <- commandArgs(TRUE)
option.list <- list(
    make_option(c("-t", "--treatTab"),       action="store", type="character", metavar="FILE",
                help="Path of treatment file for experimental peaks. This file specifies all samples
                in the project, including replicates. Required columns are ID (must be identical
                across replicates), Replicate, Batch, Barcode (e.g., W001), Plate, Treatment."),
    make_option(c("-n", "--norm"),           default="wpMedian",
                help="Type of normalization. Supported are 'none', 'pMedian' (plate median), wpMedian'
                (plate median with higher weight for neg. controls) [%default]"),
    make_option(c("-w", "--negCtlWt"),       default=20,
                help="Weight for neg. controls when --norm=wpMedian [%default]"),
    make_option(c("-m", "--minCounts"),      default=20,
                help="Minimum counts per event for filtering [%default]"),
    make_option(c("-c", "--cores"),          default=1,
                help="CPU cores to run with [%default]")
)
parser <- OptionParser(option_list=option.list,
                       usage="usage: %prog [options] OUTDIR",
                       description="Normalize raw PSI and produce dPSI and SSMD")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) != 1)
    stop(print_help(parser))


### Function definition
normalizePSI <- function(x, treat, type, negCtlWt=1) {
### Normalize raw PSI
    plates <- unique(treat$Plate)
                      
    if (type == "none") {
        norm <- as.matrix(raw[,-c(1,2)])
        rownames(norm) <- paste(raw$ID, raw$Replicate, sep=".")
    }
    
    if (type %in% c("pMedian", "pMedian")) {
        norm <- matrix(nrow=nrow(raw), ncol=ncol(raw) - 2,
                       dimnames=list(paste(raw$ID, raw$Replicate, sep="."), names(raw)[-c(1,2)])
                       )
        for (i in 1:ncol(norm)) {
            allMed <- median(x[treat$Type != "ctlPos", i + 2], na.rm=T)  
            plateMed <- sapply(plates, FUN=function(y) {
                median(c(x[treat$Type != "ctlPos" & treat$Plate == y, i + 2],
                         rep(x[treat$Type == "ctlNeg" & treat$Plate == y, i + 2], negCtlWt - 1)), na.rm=T)
            })

            no <- x[, i + 2]
            for (j in plates) {no[treat$Plate == j] <- no[treat$Plate == j] - plateMed[plates == j] + allMed}
            norm[,i] <- pmax(0, pmin(100, no))
        }
    }

    norm
}

calcSSMD <- function(x, treat, counts, opt) {  
### Calculate SSMD based on (normalized) PSI

    ## Find replicates and check
    type.neg  <- which(treat$Type == "ctlNeg")
    id.reps <- lapply(unique(treat$ID), FUN=function(y) {which(treat$ID == y)})
    names(id.reps) <- unique(treat$ID)

    nreps <- sapply(id.reps, length)
    if (any(nreps < 2)) {
        warning(length(which(nreps < 2)), " treatment ID(s) have less than 2 replicates... removing")
        id.reps <- types.reps[nreps >= 2]
    }

    ## Calculate SSMD
    ssmd <- lapply(1:ncol(x), FUN=function(y) {
        calculateSSMD(x[,y], counts[,y + 2], id.reps, type.neg, cores=opt$options$cores, ev=colnames(x)[y])
    })
    ssmd <- list(psi  = sapply(ssmd, "[[", 1),
                 dpsi = sapply(ssmd, "[[", 2),
                 ssmd = sapply(ssmd, "[[", 3)
                 )
    dimnames(ssmd$psi) <- dimnames(ssmd$dpsi) <- dimnames(ssmd$ssmd) <-
        list(names(id.reps), colnames(x))

    ## Filter treatments with counts below minimum
    for (i in 1:ncol(x)) {
        few <- which(sapply(id.reps, FUN=function(y) {any(counts[y,i + 2] < opt$options$minCounts)}))
        ssmd$psi[few,i]  <- NA
        ssmd$dpsi[few,i] <- NA
        ssmd$ssmd[few,i] <- NA

    }
    ssmd$psi  <- round(ssmd$psi, 2)
    ssmd$dpsi <- round(ssmd$dpsi, 2)
    ssmd$ssmd <- round(ssmd$ssmd, 3)
                          
    ssmd
}

calculateSSMD <- function(psi, counts, id.reps, type.neg, ev="?", cores=1) {
### Called by calcSSMD
    psi <- psi/100
    shapesNeg <- retry.posterior.shapes(psi[type.neg], mean(counts[type.neg]), ev=ev, type="neg. controls")
    psiNeg <- beta.mean(shapesNeg)
    varNeg <- beta.var(shapesNeg)

    out <- mclapply(1:length(id.reps), FUN=function(x) {
        get.ssmd(psi[id.reps[[x]]], counts[id.reps[[x]]], psiNeg, varNeg, ev=ev, type=names(id.reps)[x])
    }, mc.cores=cores)
    out <- data.frame(PSI  = 100 * sapply(out, "[[", 1),
                      SSMD = sapply(out, "[[", 2)
                      )
    out$dPSI <- out$PSI - 100 * psiNeg
    out[,c(1,3,2)]    
}

retry.posterior.shapes <- function(psi, n, ev="?", type="?", maxit=200) {
### Try running posterior.shapes() until there is a fit, but at most 'maxit' times.
    it <- 1
    out <- c(NA)
    while (it <= maxit) {
        err <- try(out <- posterior.shapes(psi, n), silent=TRUE)
        if (!is.na(out[1])) break
        if (any(grepl("optimization failed", err))) {
            it <-  it + 1
            if (it > maxit) {
                cat("[calcSSMD] EM failed after ", maxit, " tries at ", ev,
                    ", treatment ", type, "\n", sep="")
                break
            }
        } else {
            break
        }
    }
    return(out)
}

posterior.shapes <- function(psi, n) {
### For a vector of PSI, return shape parameters of a max.likelihood fitted beta distribution,
### initializing with shapes based on n.
    psi[psi < 0] <- 0
    psi[psi > 1] <- 1
    
    if (is.na(n)) {
        return(c(NA, NA))
    } else {
        psi <- psi[!is.na(psi)]
        noise <- runif(length(psi), 0, 0.00001)
        psi <- psi + ifelse(psi <= 0.5, noise, -noise)
        fit <- fitdistr(psi, "beta", start=list(shape1=1, shape2=1),
                       method="L-BFGS-B", lower=0.01, upper=n)$estimate
        round(fit, 3)
    }
}

beta.mean <- function(shapes) {
### Expected mean value of a beta distribution from shape params
  shapes[1] / sum(shapes)
}

beta.var <- function(shapes) {
### Variance of a beta distribution
    (shapes[1]*shapes[2]) / ((shapes[1]+shapes[2])^2*(shapes[1]+shapes[2]+1))
}

get.ssmd <- function(psi, counts, psiNeg, varNeg, ev="?", type="?") {
### Wrap beta.ssdm given the deeds for the controls
    shapes <- retry.posterior.shapes(psi, mean(counts), ev=ev, type=type)
    psiMean <- beta.mean(shapes)
    ssmd <- (psiMean - psiNeg) / sqrt(beta.var(shapes) + varNeg)
    c(mean = psiMean, ssmd = ssmd)
}


### Load packages
libMissing <- !require(MASS, quietly=T) && stop("Failed to load R package 'MASS'")
libMissing <- !require(parallel, quietly=T) && stop("Failed to load R package 'parallel'")


### Check input
inDir <- sub("\\/*$","",opt$args[length(opt$args)]) # directories are checked upstream

if (is.null(opt$options$treatTab))          {stop("--treatTab must be specified")}
if (!(opt$options$norm %in% c("none","pMedian","wpMedian")))  {
    stop("--norm must be one of 'none', 'pMedian', 'wpMedian'")
}

if (!dir.exists(inDir))                     {stop("Input directory not found")}
if (!dir.exists(file.path(inDir, "raw")))   {stop("Subdirectory raw/ not found")}
if (!dir.exists(file.path(inDir, "norm")))  {dir.create(file.path(inDir, "norm"))}
if (!file.exists(opt$options$treatTab))     {stop("Treatment table not found")}
if (!file.exists(file.path(inDir, "raw", "PSI.raw.tab"))) {
    stop("Raw PSI file not found at ", paste(inDir, "/raw/PSI.raw.tab", sep=""))
}
if (!file.exists(file.path(inDir, "raw", "ReadsPerEvent.tab"))) {
    stop("Raw PSI file not found at ", paste(inDir, "/raw/ReadsPerEvent.tab", sep=""))
}

### Load input
treat  <- read.delim(opt$options$treatTab)
raw    <- read.delim(file.path(inDir, "raw", "PSI.raw.tab"), check.names=F)
counts <- read.delim(file.path(inDir, "raw", "ReadsPerEvent.tab"), check.names=F)

### Normalize
norm <- normalizePSI(raw,
                     treat    = treat,
                     type     = opt$options$norm,
                     negCtlWt = switch(opt$options$norm, pMedian=1, opt$options$negCtlWt)
                     )

ssmd <- calcSSMD(norm,
                 treat  = treat,
                 counts = counts,
                 opt    = opt)


### Save tables
uniqTreat <- merge(data.frame(ID = rownames(ssmd$psi), ssmdInd = 1:nrow(ssmd$psi)),
                   data.frame(ID = treat$ID, Treatment = treat$Treatment)[!duplicated(treat$ID),],
                   by=1, all.x=T)
uniqTreat <- uniqTreat[order(uniqTreat$ssmdInd),-2]

write.table(data.frame(uniqTreat, ssmd$psi,  check.names=F), row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "norm", "PSI.norm.tab"))
write.table(data.frame(uniqTreat, ssmd$dpsi, check.names=F), row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "norm", "dPSI.norm.tab"))
write.table(data.frame(uniqTreat, ssmd$ssmd, check.names=F), row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "norm", "SSMD.norm.tab"))




##################
stop("Done")
opt <- list(options=list(), args="/media/DataExt4/SPAR-Seq_neural/full170629/seq/PSI/")
opt$options$norm = "wpMedian"
opt$options$negCtlWt = 20
opt$options$treatTab = "/media/DataExt4/SPAR-Seq_neural/full170629/seq/treatments/Treatments.fullScreen_UB170712.tab"
opt$options$cores = 4
