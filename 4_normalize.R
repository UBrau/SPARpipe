#!/usr/bin/env Rscript
### Normalize raw PSI and compute SSMD
### Part of the SPARpipe pipeline for analysis of SPAR-seq data
###
### U. Braunschweig, 08/2017


libMissing <- !require(optparse, quietly=T) && stop("Failed to load R package 'optparse'")

### Parse input
args <- commandArgs(TRUE)
option.list <- list(
    make_option(c("-t", "--treatTab"),       action="store", type="character", metavar="FILE",
                help="Path of treatment table. This file specifies all samples
                in the project, including replicates. Required columns are ID (must be identical
                across replicates), Replicate, Barcode (e.g., W001), Plate, Treatment"),
    make_option(c("-n", "--norm"),           default="pMedian",
                help="Type of normalization. Supported are 'none', 'pMedian' (plate median), wpMedian'
                (plate median with higher weight for neg. controls) [%default]"),
    make_option(c("-w", "--negCtlWt"),       default=5,
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
countsFromShapes <- function(x) {
### Translate pseudocounts in the form incl=excl to counts
    tmp <- as.numeric(unlist(strsplit(as.character(x), split="=")))
    tmp[seq(1, length(tmp) - 1, 2)] + tmp[seq(2, length(tmp), 2)]
}

normalizePSI <- function(x, treat, type, negCtlWt=1) {
### Normalize raw PSI
    plates <- unique(treat$Plate)
    if (type == "none") {
        norm <- as.matrix(raw[,-c(1,2)])
        rownames(norm) <- paste(raw$ID, raw$Replicate, sep=".")
    }
    
    if (type %in% c("pMedian", "wpMedian")) {
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
    no.reps <- nreps < 2

    if (length(type.neg) > 2) {
        stop("Unable to calculate SSMD with less than two negative controls")
    }
    if (all(no.reps)) {
        stop("No treatment has replicates, cannot calculate SSMD.")
    }
    if (any(nreps < 2)) {
        warning(length(which(nreps < 2)), " treatment ID(s) have less than 2 replicates... removing")
    }

    ## Calculate SSMD
    ssmd <- lapply(1:ncol(x), FUN=function(y) {
        calculateSSMD(x[,y], counts[,y + 2], id.reps, no.reps, type.neg,
                      cores=opt$options$cores, ev=colnames(x)[y], minCounts=opt$options$minCounts)
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

calculateSSMD <- function(psi, counts, id.reps, no.reps, type.neg, ev="?", cores=1, minCounts) {
### Called by calcSSMD
    psi <- psi/100
    shapesNeg <- retry.posterior.shapes(psi[type.neg], mean(counts[type.neg]), ev=ev, type="neg. controls",
                                        minN=minCounts)
    psiNeg <- beta.mean(shapesNeg)
    varNeg <- beta.var(shapesNeg)
    out <- mclapply(1:length(id.reps), FUN=function(x) {
    	if (no.reps[x]) {
	    return(c(mean=NA, ssmd=NA))
	} else {
            return(get.ssmd(psi[id.reps[[x]]], counts[id.reps[[x]]], psiNeg, varNeg, ev=ev, type=names(id.reps)[x],
                        minCounts=minCounts)
            )
        }
    }, mc.cores=cores)
    out <- data.frame(PSI  = 100 * sapply(out, "[[", 1),
                      SSMD = sapply(out, "[[", 2)
                      )
    out$dPSI <- out$PSI - 100 * psiNeg
    out[,c(1,3,2)]    
}

retry.posterior.shapes <- function(psi, n, ev="?", type="?", minN=1, maxit=200) {
### Try running posterior.shapes() until there is a fit, but at most 'maxit' times.
    it <- 1
    out <- c(NA)

    if (n >= minN) {
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

get.ssmd <- function(psi, counts, psiNeg, varNeg, ev="?", type="?", minCounts=1) {
### Wrap beta.ssdm given the deeds for the controls
    shapes <- retry.posterior.shapes(psi, mean(counts), ev=ev, type=type, minN=minCounts)
    psiMean <- beta.mean(shapes)
    ssmd <- (psiMean - psiNeg) / sqrt(beta.var(shapes) + varNeg)
    c(mean = psiMean, ssmd = ssmd)
}



drawPlots <- function(counts, raw, norm, ssmd.raw, ssmd.norm, opt, inDir, treat) {
### Generate plots to monitor normalization and screen performance

    ## %samples for each event with low read counts
    pdf(file.path(inDir, "norm", "LowCountEvents.pdf"), wid=3 + ncol(norm)/12, hei=5)
    plotLowCounts(counts, minCounts=opt$options$minCounts)
    dev.off()
    
    ## Sorted raw and normalized PSI
    pdf(file.path(inDir, "norm", "PSIsamples_raw.norm.pdf"), wid=9, hei=8)
    plotSamplePSI(raw, norm, treat)
    dev.off()

    ## dPSI before and after normalization
    pdf(file.path(inDir, "norm", "PSIreplicateDiff_raw.norm.pdf"), wid=6, hei=6)
    medDiff <- plotRepDiff(raw, norm, treat)
    dev.off()

    pdf(file.path(inDir, "norm", "PSIreplicateDiff_change.pdf"), wid=3 + ncol(norm)/8, hei=6)
    plotMedRepDiff(medDiff)
    dev.off()

    ## Sorted SSMD before and after normalization
    pdf(file.path(inDir, "norm", "SSMDsorted_raw.norm.pdf"), wid=6, hei=8)
    plotTreatSSMD(ssmd.raw$ssmd, ssmd.norm$ssmd, treat=treat)
    dev.off()

    ## Scatter plot of SSMD before and after normalization
    pdf(file.path(inDir, "norm", "SSMDscatter_raw.norm.pdf"), wid=6, hei=6.5)
    plotScatterSSMD(ssmd.raw$ssmd, ssmd.norm$ssmd, treat=treat)
    dev.off()
    
    ## Screening window
    if (!("ctlPos" %in% treat$Type[!duplicated(treat$ID)])) {
        cat("No positive controls found, skipping performance plot. Annotate positive and negative controls",
            "with 'ctlPos' and 'ctlNeg' in column 'Type' of the treatment table.\n")
    } else {
        pdf(file.path(inDir, "norm", "ScreenPerformance.pdf"), wid=6, hei=10)
        plotPerformance(ssmd.raw$ssmd, ssmd.norm$ssmd, treat, cores=opt$options$cores)
        dev.off()
    }
}

plotLowCounts <- function(counts, minCounts=1) {
    id.reps <- lapply(unique(treat$ID), FUN=function(y) {which(treat$ID == y)})
    names(id.reps) <- unique(treat$ID)

    lowCounts <- sapply(3:ncol(counts), FUN=function(x) {
        nLow <- length(which(sapply(id.reps, FUN=function(y) {
            any(counts[y, x] < minCounts)
        })))
    })

    par(mar=c(4.3, 4, 3.5, 2))
    xpos <- barplot(rep(NA, length(lowCounts)), yaxt="n", ylim=c(0,130), ylab="% of samples",
            main=paste("Treatments with less than ", minCounts, " reads in at least one replicate", sep=""),
            xlab=paste(length(which(lowCounts/length(id.reps) < 0.1)), "/", length(lowCounts), 
                       " events with more than 90% treatments above 20 reads", sep="")
            )
    abline(h=seq(0,100,10), col="grey90")
    axis(2, at=seq(0, 100, 20))
    barplot(100 * lowCounts/length(id.reps), ylim=c(0,100), add=T, yaxt="n", )
    par(xpd=NA)
    text(xpos, 100 * lowCounts/length(id.reps), srt=90, adj=c(-0.08, 0.5), labels=names(counts)[-c(1,2)], cex=0.7)
    par(xpd=FALSE)
}

plotSamplePSI <- function(raw, norm, treat) {
### Plot all raw and normalized PSI in a row, one page per event
    cols <- rep("black", nrow(treat))
    cols[treat$Type == "exp"] <- "grey70"
    cols[treat$Type == "ctlNeg"] <- "dodgerblue"
    cols[treat$Type == "ctlPos"] <- "brown1"
    pch <- ifelse(grepl("ctl", treat$Type), 20, 1)
    plateConsec <- length(which(diff(as.integer(treat$Plate)) != 0)) == (length(unique(treat$Plate)) - 1)  # plates are all together
    plateLen    <- table(treat$Plate)

    par(mfrow=c(2,1), mar=c(4.2,4,3,2))
    for (i in 1:ncol(norm)) {
        suppressWarnings(ylim <- range(c(raw[,i + 2], norm[,i]), na.rm=T))
        if (any(is.infinite(ylim))) {ylim <- c(0,100)}
        .plotSamplePSI.event(raw[,i + 2], event=colnames(norm)[i], treat=treat, cols=cols, pch=pch, ylim=ylim,
                          plateConsec=plateConsec, plateLen=plateLen)
        title(names(raw)[i + 2], adj=0, col.main="black", cex.main=1.5)
        title("Raw PSI", adj=1, col.main="grey70")

        .plotSamplePSI.event(norm[,i], event=colnames(norm)[i], treat=treat, cols=cols, pch=pch, ylim=ylim,
                          plateConsec=plateConsec, plateLen=plateLen)
        title("Normalized PSI", adj=1, col.main="grey70")
    }
}

.plotSamplePSI.event <- function(psi, event, treat, cols, pch, plateConsec, plateLen, ylim=NULL, ylab="PSI") {    
    plot(psi, col=cols, pch=ifelse(grepl("ctl", treat$Type), NA, 1), xlab="Sample", ylab=ylab, ylim=ylim)
    points(psi, col=cols, pch=ifelse(grepl("ctl", treat$Type), 20, NA))
    if (plateConsec) {  # Plate medians
        plateMed <- sapply(unique(treat$Plate), FUN=function(x) {median(psi[treat$Plate == x], na.rm=T)})
        segments(x0=cumsum(plateLen) - (plateLen - 1), x1=cumsum(plateLen), y0=plateMed, lwd=2, col="grey30")
    }
    if (length(which(diff(as.integer(treat$Replicate)) != 0)) == (length(unique(treat$Replicate)) - 1) &
        length(unique(treat$Replicate)) > 1) {  # Replicate breaks
        abline(v=(cumsum(table(treat$Replicate)) - (table(treat$Replicate) - 0.5))[-1], lty=3, col="grey80")
    }
    
}


plotRepDiff <- function(raw, norm, treat, colR="saddlebrown", colN="yellowgreen") {
### Plot raw and norm PSI difference between replicates for all events

    id.reps <- lapply(unique(treat$ID), FUN=function(y) {which(treat$ID == y)})
    names(id.reps) <- unique(treat$ID)
    pDiff <- rep(NA, ncol(norm))
    names(pDiff) <- colnames(norm)

    for (i in 1:ncol(norm)) {
        pDiff.raw  <- .plotRepDiff.range(raw[,i + 2], id.reps)
        pDiff.norm <- .plotRepDiff.range(norm[,i], id.reps)

        if (length(which(!is.na(pDiff.raw))) > 2) {
            showR <- TRUE
            densR <- density(pDiff.raw, cut=0, na.rm=T)
            medR <- median(pDiff.raw, na.rm=T)
            legR <- paste("Raw (median=", round(medR, 2), ")", sep="")
        } else {
            showR <- FALSE
            legR <- "Raw"
        }
        if (length(which(!is.na(pDiff.norm))) > 2) {
            showN <- TRUE
            densN <- density(pDiff.norm, cut=0, na.rm=T)
            medN <- median(pDiff.norm, na.rm=T)
            legN <- paste("Normalized (median=", round(medN, 2), ")", sep="")
            pDiff[i] <- medN - medR
        } else {
            showN <- FALSE
            densN <- list(y=NA)
            legN <- "Normalized"
        }
        
            
        if (showR) {
            plot(densR, lwd=3, col="saddlebrown", xlab="PSI", main="", ylim=c(0, max(c(densR$y, densN$y), na.rm=T)))
            segments(x0=min(densR$x), y0=0, y1=densR$y[1], col=colR, lwd=3)
            segments(x0=max(densR$x), y0=0, y1=densR$y[length(densR$y)], col=colR, lwd=3)
            abline(v=medR, lty=2, col=colR)
            if (showN) {
                lines(densN <- density(pDiff.norm, cut=0, na.rm=T), lwd=3, col=colN)
                segments(x0=min(densN$x), y0=0, y1=densN$y[1], col=colN, lwd=3)
                segments(x0=max(densN$x), y0=0, y1=densN$y[length(densN$y)], col=colN, lwd=3)
                abline(v=medN, lty=2, col=colN)
            } else {
                legN <- "Normalized"
            }
        } else {
            plot(1, 1, type="n", xlim=c(0,100), ylim=c(0,1), xlab="PSI", ylab="Density")
            text(50, 0.5, "Not enough data")
        }
        title(names(raw)[i + 2], adj=0, col.main="black", cex.main=1.5)
        title("Replicate PSI difference", adj=1, col.main="grey70")
        legend("topright", lwd=3, col=c(colR, colN), legend=c(legR, legN), bty="n")
    }

    pDiff
}

.plotRepDiff.range <- function(psi, id.reps) {
    sapply(id.reps, FUN=function(x) {
        out <- suppressWarnings(range(psi[x], na.rm=T))
        if (any(is.infinite(out))) {out <- c(NA, NA)}
        out[2] - out[1]
    })
}

plotMedRepDiff <- function(medDiff, cols=c("saddlebrown","yellowgreen")) {
### Plot the change in median PSI difference
    medDiff <- sort(medDiff)
    par(mar=c(10, 5, 3.5, 3))
    xpos <- barplot(medDiff, names.arg=NA, col=ifelse(medDiff <= 0, cols[2], cols[1]),
                    ylab="Change in median replicate PSI difference\n(norm-raw)")
    par(xpd=T)
    text(xpos, ifelse(medDiff <= 0, medDiff, 0), srt=90, adj=c(1.1, 0.5), labels=names(medDiff))
    par(xpd=NA)
    title(main="Change in replicate PSI difference during normalization", cex.main=1.2)    
}

plotTreatSSMD <- function(raw, norm, treat) {
### Plot all raw and normalized PSI in a row, one page per event
    treat1 <- treat[!duplicated(treat$ID),]
    cols <- rep("black", nrow(treat1))
    cols[treat1$Type == "exp"] <- "grey70"
    cols[treat1$Type == "ctlNeg"] <- "dodgerblue"
    cols[treat1$Type == "ctlPos"] <- "brown1"
    pch <- ifelse(grepl("ctl", treat1$Type), 20, 1)
    plateConsec <- length(which(diff(as.integer(treat1$Plate)) != 0)) == (length(unique(treat1$Plate)) - 1)  # plates are all together
    plateLen    <- table(treat1$Plate)

    par(mfrow=c(2,1), mar=c(4.2,4,3,2))
    for (i in 1:ncol(norm)) {
        suppressWarnings(ylim <- range(c(raw[,i], norm[,i]), na.rm=T))
        if (any(is.infinite(ylim))) {ylim <- c(-3,3)}
        .plotTreatSSMD.event(raw[,i], event=colnames(norm)[i], treat=treat1, cols=cols, pch=pch, ylim=ylim,
                          plateConsec=plateConsec, plateLen=plateLen, ylab="SSMD")
        title(colnames(raw)[i], adj=0, col.main="black", cex.main=1.5)
        title("Raw SSMD", adj=1, col.main="grey70")

        .plotTreatSSMD.event(norm[,i], event=colnames(norm)[i], treat=treat1, cols=cols, pch=pch, ylim=ylim,
                          plateConsec=plateConsec, plateLen=plateLen, ylab="SSMD")
        title("Normalized SSMD", adj=1, col.main="grey70")
    }
}

.plotTreatSSMD.event <- function(psi, event, treat, cols, pch, plateConsec, plateLen, ylim=NULL, ylab="PSI") {    
    plot(psi, col=cols, pch=ifelse(grepl("ctl", treat$Type), NA, 1), xlab="Sample", ylab=ylab, ylim=ylim)
    abline(h=c(-3,0,3), lty=c(2,1,2), col="grey30")
    points(psi, col=cols, pch=ifelse(grepl("ctl", treat$Type), 20, NA))
    if (plateConsec) {  # Plate medians
        plateMed <- sapply(unique(treat$Plate), FUN=function(x) {median(psi[treat$Plate == x], na.rm=T)})
        segments(x0=cumsum(plateLen) - (plateLen - 1), x1=cumsum(plateLen), y0=plateMed, lwd=2, col="grey30")
    }
}

plotScatterSSMD <- function(raw, norm, treat) {
    treat1 <- treat[!duplicated(treat$ID),]
    cols <- rep("black", nrow(treat1))
    cols[treat1$Type == "exp"] <- "grey70"
    cols[treat1$Type == "ctlNeg"] <- "dodgerblue"
    cols[treat1$Type == "ctlPos"] <- "brown1"
    pch <- ifelse(grepl("ctl", treat1$Type), 20, 1)

    negRangeCol <- tCols("dodgerblue", transp=0.08)
    posRangeCol <- tCols("brown1", transp=0.08)

    for (i in 1:ncol(norm)) {
        if (length(which(!is.na(norm[,i]))) > 0) {
            lims <- range(c(raw[,i], norm[,i]), na.rm=T)
            suppressWarnings(rangeNegR <- range(raw[treat1$Type  == "ctlNeg",i], na.rm=T))
            suppressWarnings(rangeNegN <- range(norm[treat1$Type == "ctlNeg",i], na.rm=T))
            suppressWarnings(rangePosR <- range(raw[treat1$Type  == "ctlPos",i], na.rm=T))
            suppressWarnings(rangePosN <- range(norm[treat1$Type == "ctlPos",i], na.rm=T))
            
            plot(raw[,i], norm[,i], xlim=lims, ylim=lims, xlab="SSMD (raw)", ylab="SSMD (norm)", type="n")
            add <- (lims[2] - lims[1])*0.1

            if (!any(is.na(c(rangeNegR, rangeNegN))) & !any(is.na(c(rangeNegR, rangeNegN)))) {
                polygon(x=c(lims[1] - add, rangeNegR[2], rangeNegR[2], lims[1] - add),
                        y=rep(rangeNegN, each=2),
                        col=negRangeCol, border=NA)
                polygon(x=c(rangeNegR[1], rep(rangeNegR[2], 2), rangeNegR[1]),
                        y=c(rep(lims[1] - add, 2), rep(rangeNegN[1], 2)),
                        col=negRangeCol, border=NA)
            }
            if (!any(is.na(c(rangePosR, rangePosN))) & !any(is.na(c(rangePosR, rangePosN)))) {
                polygon(x=c(lims[1] - add, rangePosR[2], rangePosR[2], lims[1] - add),
                        y=rep(rangePosN, each=2),
                        col=posRangeCol, border=NA)
                polygon(x=c(rangePosR[1], rep(rangePosR[2], 2), rangePosR[1]),
                        y=c(rep(lims[1] - add, 2), rep(rangePosN[1], 2)),
                        col=posRangeCol, border=NA)
            }
            
            points(raw[,i], norm[,i], pch=ifelse(pch == 1, 1, NA), col=cols)
            points(raw[,i], norm[,i], pch=ifelse(pch == 20, 20, NA), col=cols)
            abline(h=c(-3,0,3), col="grey30", lty=c(3,1,3))
            abline(v=c(-3,0,3), col="grey30", lty=c(3,1,3))
            abline(a=0, b=1, lty=2, col="grey30")
        } else {
            plot(1, 1, type="n", xlim=c(-3,3), ylim=c(-3,3), xlab="SSMD (raw)", ylab="SSMD (norm)")
            text(0, 0, "Not enough data")
        }
        title(main=paste("SSMD in", colnames(norm)[i]))
    }
}

tCols <- function(x, transp=0.1) {
    tColors <- t(col2rgb(x))
    rgb(tColors, alpha=round(255 * transp), maxColorValue=255)
}

plotPerformance <- function(raw, norm, treat, cores) {
    treat1 <- treat[!duplicated(treat$ID),]
    groups <- list("Pos. controls" = which(treat1$Type == "ctlPos"),
                   "Neg. controls" = which(treat1$Type == "ctlNeg"),
                   "Experimental"  = which(treat1$Type == "exp")
                   )
    cols <- c("brown1","dodgerblue","black")
    suppressWarnings(
        xlim <- c(0, max(c(max(abs(c(as.numeric(raw[unlist(groups[2]),]), as.numeric(norm[unlist(groups[1]),])))),
                           quantile(apply(abs(cbind(raw, norm))[groups[[1]],], MAR=1, max, na.rm=T), 0.50, na.rm=T),
                           quantile(apply(abs(cbind(raw, norm))[groups[[3]],], MAR=1, max, na.rm=T), 0.99, na.rm=T)
                           ), na.rm=T)
                  )
    )

    par(mfrow=c(2,1), mar=c(4,4.2,3,1))
    .plotPerformance.one(raw, cols=cols, groups=groups, xlim=xlim, main="Raw data", cores=cores)
    title(main="Screen performance", adj=0, col.main="black", cex.main=1.5, line=1.3)
    .plotPerformance.one(norm, cols=cols, groups=groups, xlim=xlim, main="Normalized data", cores=cores)
}

.plotPerformance.one <- function(x, cols, groups, xlim=c(0, 30), minEvents=0.25, main="", cores=1) {
### Generate the performance plot for one SSMD table
    s.thres <- seq(xlim[1], xlim[2]*1.05, 0.01)
    black <- which(apply(x, MAR=1, FUN=function(y) {length(which(is.na(y))) < minEvents * ncol(x)}))
    val <- sapply(groups, FUN=function(y) {
        y <- setdiff(y, black)
        unlist(mclapply(s.thres, .plotPerformance.count, x=x, y=y, mc.cores=cores))
    })

    plot(s.thres, val[,1], type="n", xlim=xlim, ylim=c(0, 1),
         xlab="|SSMD| threshold for hit", ylab="Fraction of hits")
    title(main=main, col.main="grey70", adj=1)
    abline(h=c(0,1))
    for (i in ncol(val):1) {lines(s.thres, val[,i], lwd=3, col=cols[i])}
    accrange <- s.thres[c(min(which(val[,2] == 0)), max(which(val[,1] == 1)))]
    abline(v=accrange, col="black", lty=2)
    legend("topright", colnames(val), text.col=cols, bty="n")

    arrwid <- (xlim[2] - xlim[1]) / 50
    segments(x0=accrange[1] + 0.5*arrwid,  x1=accrange[2] - 0.5*arrwid, y0=1.02, lwd=2)
    polygon(x=c(accrange[1], rep(accrange[1]+arrwid, 2)), y=c(1.02, 1.01, 1.03), col="black")
    polygon(x=c(rep(accrange[2] - arrwid, 2), accrange[2]), y=c(1.03, 1.01, 1.02), col="black")
    par(xpd=NA)
    text(accrange, 1.06, labels=accrange, adj=c(0.5, 0))
    par(xpd=FALSE)
}

.plotPerformance.count <- function(x, y, z) {
    if (length(y) > 1) {	      
        return(sum(apply(x[y,], MAR=1, FUN=function(a) {any(abs(a) >= z)}), na.rm=T) / length(y))
    }
    if (length(y) == 1) {
        return(sum(any(abs(x[y,]) >= z)))
    }
    if (length(y) == 0) {
        return(NA)
    }
}
#### Main part ####

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
if (!file.exists(file.path(inDir, "raw", "InclExclCounts.tab"))) {
    if (file.exists(file.path(inDir, "raw", "ReadsPerEvent.tab"))) {
        counts <- read.delim(file.path(inDir, "raw", "ReadsPerEvent.tab"))
        for (i in 3:ncol(counts)) {counts[is.na(counts[,i]), i] <- 0}
    } else {
        stop("Reads per event-file not found at ", paste(inDir, "/raw/ReadsPerEvent.tab", sep=""))
    }
} else {
    cat("Reading event read counts from InclExclCounts.tab\n")
    counts <- read.delim(file.path(inDir, "raw", "InclExclCounts.tab"), check.names=F)
    evNames <- names(counts)[-c(1,2)]
    suppressWarnings(
        counts <- data.frame(counts[,1:2], sapply(3:ncol(counts), FUN=function(x) {countsFromShapes(counts[,x])}))
    )
    for (i in 3:ncol(counts)) {counts[is.na(counts[,i]), i] <- 0}
    names(counts)[-c(1,2)] <- evNames
}

### Load input
treat  <- read.delim(opt$options$treatTab)
raw    <- read.delim(file.path(inDir, "raw", "PSI.raw.tab"), check.names=F)


### Normalize
cat("Normalizing...\n")
norm <- normalizePSI(raw,
                     treat    = treat,
                     type     = opt$options$norm,
                     negCtlWt = switch(opt$options$norm, pMedian=1, opt$options$negCtlWt)
                     )

ssmd.raw  <- calcSSMD(raw[-c(1,2)],
                      treat  = treat,
                      counts = counts,
                      opt    = opt)

ssmd.norm <- calcSSMD(norm,
                      treat  = treat,
                      counts = counts,
                      opt    = opt)


## Save tables
cat("Saving tables...\n")
uniqTreat <- merge(data.frame(ID = rownames(ssmd.norm$psi), ssmdInd = 1:nrow(ssmd.norm$psi)),
                   data.frame(ID = treat$ID, Treatment = treat$Treatment)[!duplicated(treat$ID),],
                   by=1, all.x=T)
uniqTreat <- uniqTreat[order(uniqTreat$ssmdInd),-2]

write.table(data.frame(uniqTreat, ssmd.norm$psi,  check.names=F), row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "norm", "PSI.norm.tab"))
write.table(data.frame(uniqTreat, ssmd.norm$dpsi, check.names=F), row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "norm", "dPSI.norm.tab"))
write.table(data.frame(uniqTreat, ssmd.norm$ssmd, check.names=F), row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "norm", "SSMD.norm.tab"))


### Generate plots
cat("Generating plots...\n")
drawPlots(counts, raw, norm, ssmd.raw, ssmd.norm, opt, inDir, treat)

cat("Done\n\n")
