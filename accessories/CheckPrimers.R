#!/usr/bin/env Rscript
### Given a table of forward and reverse event primers, check which ones would have a 3'-end
### overlap of at least 5 nt with which others.
###
### U. Braunschweig 2018-2021


cArgs <- commandArgs(TRUE)

## Check dependencies
libMissing <- !require("optparse", quietly=T)
if (libMissing) {
    install.packages("optparse", repos="https://cloud.r-project.org")
    libMissing <- !require("optparse", quietly=T)
    if (libMissing) {stop("Failed to load R package 'optparse'")}
}

## Capture input
opt.list <- list(
    make_option(c("-p", "--primerFile"),  action="store", type="character",
                help="Path of CSV file. Must contain columns event, seqF, seqR (can include adapters, 5'-3')"),
    make_option(c("-m", "--minOl"),       action="store", default=5,
                help="Minimum 3' end overlap to be considered critical [default: %default]"),
    make_option("--outName",              action="store", type="character",
                help="Name for output [default: will be take from input]")
)

opt <- parse_args(OptionParser(option_list=opt.list), args=cArgs)

### Function defs
revcomp <- function(x) {
    sapply(strsplit(chartr("ACGT", "TGCA", x), split=NULL),
           FUN=function(y) {paste(rev(y), collapse="")})
}


### Check and load input
if (is.null(opt$primerFile))      {stop("PRIMERFILE must be provided")}
if (!file.exists(opt$primerFile)) {stop("Could not find PRIMERFILE, file: ", opt$primerFile)}

prim <- read.csv(opt$primerFile, as.is=T, comment.char="#")
prim <- prim[!grepl("^#", prim[,1]),]

missCol <- data.frame(col   = c("event", "seqF", "seqR"),
                      found = sapply(c("event", "seqF", "seqR"), FUN=function(x) {x %in% names(prim)})
                      )
if (any(!missCol$found))          {stop("Column(s) ", paste(missCol$col[!missCol$found], collapse=", "),
                                        " not found")}
prim$seqF <- toupper(gsub(" ", "", prim$seqF))
prim$seqR <- toupper(gsub(" ", "", prim$seqR))


### Do stuff
seq <- data.frame(
    primer = paste(rep(prim$event, 2), rep(c("Fwd","Rev"), each=nrow(prim)), sep="."),
    event  = rep(prim$event, 2),
    seq    = c(prim$seqF, prim$seqR),
    stringsAsFactors=F
)

olwhich <- list()
olcount <- matrix(0, nrow=nrow(seq), ncol=min(nchar(seq$seq)), dimnames=list(seq$primer, c()))
found <- TRUE
k <- opt$minOl

while (found) {
    kmersF <- sub(paste(".*(.{", k, "})$", sep=""), "\\1", seq$seq)
    kmersR <- revcomp(kmersF)
    
    olcount[,k] <- sapply(kmersF, FUN=function(x) {length(which(kmersR == x))})
    olwhich.k <- lapply(kmersF, FUN=function(x) {seq$primer[kmersR == x]})
    if (all(sapply(olwhich.k, length) == 0))  {
        found <- FALSE
        cat("Maximum match length:", k - 1, "\n")
        break
    }
    olwhich <- c(olwhich, list(olwhich.k))
    names(olwhich)[length(olwhich)] <- k
    
    k  <- k + 1
    if (k > ncol(olcount)) found <- FALSE
}


matches <- lapply(1:nrow(seq), FUN=function(x) {
    sort(unique(unlist(lapply(olwhich, FUN="[[", x))))
})
seq$Nmatches <- sapply(matches, length)
seq$maxMatch <- apply(olcount, MAR=1, FUN=function(x) {
    out <- suppressWarnings(max(which(x > 0)))
    ifelse(is.infinite(out), NA, out)
})
seq$whichMatch    <- sapply(matches, paste, collapse=",")
seq$whichMaxMatch <- sapply(1:nrow(seq), FUN=function(x) {
    if (is.na(seq$maxMatch[x])) {
        return("")
    } else {
        return(paste(sort(olwhich[[seq$maxMatch[x] - opt$minOl + 1]][[x]]), collapse=","))
    }
})


### Output 

out <- seq[order(-seq$Nmatches, -seq$maxMatch, seq$primer),]
if (is.null(opt$outName)) {
    outName <- sub("\\.[^.]*$", "_primerHybrids.csv", opt$primerFile)
} else {
    outName <- opt$outName
}
write.csv(out, file=outName, row.names=F)
cat("Wrote", outName, "\n")

cat(length(which(seq$Nmatches > 0)), "primer(s) from",
    length(unique(seq$event[which(seq$Nmatches > 0)])), "pair(s) have matches of", opt$minOl, "or more bp\n")

