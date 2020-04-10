#!/usr/bin/env Rscript
### Generate junction FASTA files for SPAR-Seq screen based on a CSV with coordinates and sequences
### of every exon element in the amplicon.
### Columns are: gene, event, structure, chrom, strand, C1.start, C1.end, C2.start, C2.end, and
### start and end for each optional element (E1:En).
### The string in the 'structure' column defines which elements are cassette exons, alternative ss,
### introns and 'constitutive' introns and for which elements PSI is requested.
###
### U. Braunschweig 2017-2020
### Changes: - Sanity checks for input files


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
    make_option(c("-e", "--eventFile"),   action="store", type="character",
                help="Path of CSV file. Must contain columns gene, event, structure, chrom, strand
                and start, end and .start/.end/.seq for C1, C2, and alternative segments E1...En
                (filled starting from the lowest) if applicable.
                The string in the 'structure' column defines which elements are cassette exons,
                alternative ss, introns and 'constitutive' introns and for which elements PSI is requested.
                See https://github.com/UBrau/SPARpipe."),
    make_option(c("-p", "--primerFile"),  action="store", type="character",
                help="Path of CSV file. Must contain columns event, seqF, seqR (including adapters, 5'-3')"),
    make_option("--adaptF",                action="store", default="ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                help="Common adapter sequence in fwd primers, 5'-3' [default: %default]"),
    make_option("--adaptR",                action="store", default="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                help="Common adapter sequence in rev primers, 5'-3' [default: %default]"),
    make_option(c("-l", "--readLength"),  action="store", default=150,
                help="Read length used in sequencing run [default: %default]"),
    make_option("--trimR1.5",             action="store", default=0,
                help="Nucleotides to trim off the 5'-end of fwd reads [default: %default]"),
    make_option("--trimR1.3",             action="store", default=0,
                help="Nucleotides to trim off the 3'-end of fwd reads [default: %default]"),
    make_option("--trimR2.5",             action="store", default=0,
                help="Nucleotides to trim off the 5'-end of rev reads [default: %default]"),
    make_option("--trimR2.3",             action="store", default=0,
                help="Nucleotides to trim off the 3'-end of rev reads [default: %default]"),
    make_option("--minCnt",              action="store", default=4,
                help="Min nt of fwd/rev reads that must overlap with the C1/C2 [default: %default]"),
    make_option("--eventOut",            action="store", default="EventJunctions.tab",
                help="Name for event table output [default: %default]"),
    make_option("--juncOut",             action="store", default="Junctions",
                help="Base name for junction library FASTA output [default: %default]")
)

opt <- parse_args(OptionParser(
    usage = "Generate junction and amplicon FASTA files for SPAR-Seq screen based on a CSV file
       with coordinates and sequences of every exon element in the amplicon and a primer CSV file.
       Junctions will also be checked for minimum Hamming distance to other junctions in the same gene
       as well as all other junctions, separately for fwd and rev junctions. 
       Also produced are BED files and an event table that will be required in step 3_combine.pl.

           MakeJunctionsFASTA.R -e EVENTFILE -p PRIMERFILE [OPTIONS]",
    option_list=opt.list), args=cArgs)


### Function definitions
minDistJunc <- function(juncFA) {
### Determine minimal distance between each sequence and the other isoforms of
### the same event as well as all other sequences.
    require(Biostrings)
    
    faNam <- sub("> *", "", juncFA[seq(1, length(juncFA) - 1, 2)])
    faSeq <- juncFA[seq(2, length(juncFA), 2)]
    distMat <- as.matrix(stringDist(faSeq, method="hamming"))
    sameEvent <- as.numeric(as.factor(sub("([^_]+)_[fwrv]+[0-9]+_.+", "\\1", faNam)))
    i = 1

    distEvent <- sapply(1:length(faNam), FUN=function(x) {
        eventRows <- which(sameEvent == sameEvent[x])
        if (length(eventRows) == 1) {
            return(NA)
        } else {
            return(min(c(distMat[x, setdiff(eventRows, x)], distMat[setdiff(eventRows, x), x])))
        }        
    })
    distAll <- sapply(1:length(faNam), FUN=function(x) {
        min(c(distMat[x, -x], distMat[-x, x]))
    })

    out <- data.frame(junction = faNam, distEvent, distAll)
    out[order(out$junction),]
}



### Check and load input
if (is.null(opt$eventFile))       {stop("EVENTFILE must be provided")}
if (is.null(opt$primerFile))      {stop("PRIMERFILE must be provided")}
if (!file.exists(opt$eventFile))  {stop("Could not find EVENTFILE, file: ", opt$eventFile)}
if (!file.exists(opt$primerFile)) {stop("Could not find PRIMERFILE, file: ", opt$primerFile)}


## Check eventFile
events <- read.csv(opt$eventFile, as.is=T)
minEventCols <- c("gene","event","structure","chrom","strand","C1.start","C1.end","C2.start","C2.end","C1.seq","C2.seq")
if (!all(minEventCols %in% names(events))) {
    stop("eventFile must contain at least columns ", paste(minEventCols, collapse=", "))
}
nameFail <- grepl("_", events$event)
if (any(nameFail)) {
    stop("Event(s) ", paste(events$event[nameFail], collapse=", "), " contain(s) '_'.")
}
duploEvents <- unique(events$event[duplicated(events$event)])
if (length(duploEvents) > 0) {
    stop("Duplicated event designation(s): ", paste(duploEvents, collapse=", "))
}

for (i in grep("\\.seq", names(events))) {events[,i] <- toupper(events[,i])}


## Check primerFile
primers <- read.csv(opt$primerFile, as.is=T)
minPrimerCols <- c("event","seqF","seqR")
if (!all(minPrimerCols %in% names(primers))) {
    stop("primerFile must contain at least columns ", paste(minPrimerCols, collapse=", "))
}

if (!all(events$event == primers$event)) {stop("Columns 'event' in events and primers files do not match")}

libMissing <- !require("Biostrings", quietly=T)
if (libMissing) {stop("Package 'Biostrings' could not be loaded")}

primers$seqF <- toupper(gsub(" ","",primers$seqF))
primers$seqR <- toupper(gsub(" ","",primers$seqR))
primers$specF <- sapply(as.character(primers$seqF), FUN=function(x) {strsplit(x, split=opt$adaptF)[[1]][2]})
primers$specR <- sapply(as.character(primers$seqR), FUN=function(x) {strsplit(x, split=opt$adaptR)[[1]][2]})
primers$specR <- as.character(reverseComplement(DNAStringSet(primers$specR)))


## Check if all primers match their C1 or C2 sequences, and where they map
combi <- list(event = events$event)
combi$match.F = lapply(1:nrow(events), FUN=function(x) {
    seqs <- c(events$C1.seq[x], as.character(events[x,grep("E[0-9]+\\.seq", names(events))]))
    seqs <- seqs[!is.na(seqs)]
    tmp <- sapply(1:length(seqs), FUN=function(y) {paste(seqs[1:y], collapse="")})
    if (length(seqs) > 1) {
        tmp <- c(tmp, paste(seqs[-1], collapse=""))
        if (grepl(primers$specF[x], tmp[length(tmp)])) {stop("Fwd primer does not overlap with C1 in ", events$event[x])}
    }
    matched <- grep(primers$specF[x], tmp)
    if (length(matched) == 0) {stop("Fwd primer does not match in ", events$event[x])}
    out <- 1:min(matched)
    if (length(out) > 1) {warning("Fwd primer in ", events$event[x], " maps across first ", length(out),
                                  " elements. These will occur in all fwd junctions.")}
    out
})
combi$match.R = lapply(1:nrow(events), FUN=function(x) {
    seqs <- c(as.character(events[x,grep("E[0-9]+\\.seq", names(events))]), events$C2.seq[x])
    seqs <- seqs[!is.na(seqs)]
    tmp <- sapply(length(seqs):1, FUN=function(y) {paste(seqs[y:length(seqs)], collapse="")})
    if (length(seqs) > 1) {
        tmp <- c(tmp, paste(seqs[-length(seqs)], collapse=""))
        if (grepl(primers$specR[x], tmp[length(tmp)])) {stop("Rev primer does not overlap with C2 in ", events$event[x])}
    }
    matched <- grep(primers$specR[x], tmp)
    if (length(matched) == 0) {stop("Rev primer does not match in ", events$event[x])}
    out <- (length(seqs) + 2 -min(matched)):(length(seqs) + 1)
    if (length(out) > 1) {warning("Rev primer in ", events$event[x], " maps across last ", length(out),
                                  " elements. These will occur in all rev junctions.")}
    out
})


## Make object containing all possible sequence elements
## Generate the snippet of C1 after and of C2 before and including the primer binding site.
## Since some primers map across the junctions, first obtain the snippet distal of the primer to subtract.

combi$C1 <- sapply(1:nrow(events), FUN=function(x) {
    tmp  <- setdiff(as.character(events[x, grep("E[0-9]+\\.seq", names(events))]), NA)
    snip <- strsplit(paste(c(events$C1.seq[x], tmp), collapse=""), split=primers$specF[x])[[1]][1]
    if (is.na(snip)) {snip <- ""}
    substr(events$C1.seq[x], start = nchar(snip) + 1, stop = nchar(events$C1.seq[x]))
})

combi$C2 <- sapply(1:nrow(events), FUN=function(x) {
    tmp  <- setdiff(as.character(events[x, grep("E[0-9]+\\.seq", names(events))]), NA)
    snip <- strsplit(paste(c(tmp, events$C2.seq[x]), collapse=""), split=primers$specR[x])[[1]][2]
    if (is.na(snip)) {snip <- ""}
    substr(events$C2.seq[x], start = 1, stop = nchar(events$C2.seq[x]) - nchar(snip))
})

combi$E <- lapply(1:nrow(events), FUN=function(x) {
    out <- as.character(events[x, grep("E[0-9]+\\.seq", names(events))])
    out[!is.na(out)]
})


## Check if there are still at least minCnt nt left of C1 and C2 after 5'-trimming
trim5len.C1 <- sapply(combi$C1, nchar) - opt$trimR1.5
trim5len.C2 <- sapply(combi$C2, nchar) - opt$trimR2.5

if (all(trim5len.C1 >= opt$minCnt)) {
    cat("Requested R1 5' trimming OK. You could trim a max of", min(sapply(combi$C1, nchar)) - opt$minCnt,
        "nt to maintain C1 overlap.\n")
} else {
    stop("Trimming of ", opt$trimR1.5, " nt at R1 results in too short C1 for ",
         paste(events$event[trim5len.C1 < opt$minCnt], collapse=", "), ". Max is ",
         min(sapply(combi$C1, nchar)) - opt$minCnt, ".")
}
if (all(trim5len.C2 >= opt$minCnt)) {
    cat("Requested R2 5' trimming OK. You could trim a max of", min(sapply(combi$C2, nchar)) - opt$minCnt,
        "nt to maintain C2 overlap.\n")
} else {
    stop("Trimming of ", opt$trimR2.5, " nt at R2 results in too short C2 for ",
         paste(events$event[trim5len.C2 < opt$minCnt], collapse=", "), ". Max is ",
         min(sapply(combi$C2, nchar)) - opt$minCnt, ".")
}


## Prepare things in common for fwd and rev libraries
readEndUp <- opt$readLength - opt$trimR1.3
readEndDn <- opt$readLength - opt$trimR2.3


## Generate fwd junctions
forward <- lapply(1:nrow(events), FUN=function(x) {
    elem  <- strsplit(events$structure[x], split="[-:]")[[1]]
    eName <- gsub( "\\[.\\]", "", elem)
    eType <- {
        out <- ifelse(grepl("\\[", elem), "", "A")
        out[grepl("\\[3\\]", elem)] <- "3"
        out[grepl("\\[5\\]", elem)] <- "5"
        out[grepl("\\[i\\]", elem)] <- "I"
        out[grepl("\\[c\\]", elem)] <- "C"
        out[c(1, length(out))]      <- "C"
        out
    }
    nE <- length(combi$E[[x]])
  
    Eind <- matrix(ncol=nE + 2, nrow=2^nE)
    Eind[,c(1, ncol(Eind))] <- TRUE  
    if (length(combi$E[[x]]) > 0) {
        for (i in 1:nE) {
            Eind[,i + 1] <- rep(rep(c(FALSE, TRUE), each=2^(i-1)), times=2^nE / (2^i))
        }
    }

    lengthsAll <- c(nchar(combi$C1[x]), nchar(combi$E[[x]]), nchar(combi$C2[x]))
    lengths <- lapply(1:nrow(Eind), FUN=function(y) {cumsum(lengthsAll[Eind[y,]])})
    select <- sapply(lengths, FUN=function(y) {  # how many elements needed?
        suppressWarnings(out <- min(which(y >= readEndUp)))
        if (is.infinite(out)) {out <- length(y)}
        out
    })
    
    jNames <- sapply(1:nrow(Eind), FUN=function(y) {  # junction names using only the necessary elements
        paste(eName[which(Eind[y,])[1:select[y]]], collapse="")})
    seqsAll <- c(combi$C1[x], combi$E[[x]],combi$C2[x])
    seqs <- sapply(1:nrow(Eind), FUN=function(y) {   # seqs of the necessary elements
        paste(seqsAll[which(Eind[y,])[1:select[y]]], collapse="")})
      
    ## Remove 'forbidden' combinations due to element type and primer overlap
    forbType <- rep(FALSE, nrow(Eind))
    link5 <- which(eType %in% c("5","I"))
    if (length(link5) > 0) {
        forbType <- forbType | apply(sapply(link5, FUN=function(y) {Eind[,y] & !Eind[,y - 1]}), MAR=1, any)
    }
    link3 <- which(eType %in% c("3","I"))
    if (length(link3) > 0) {
        forbType <- forbType | apply(sapply(link3, FUN=function(y) {Eind[,y] & !Eind[,y + 1]}), MAR=1, any)
    }
    const <- which(eType == "C")
    if (length(const) > 2) {
        forbType <- forbType | !apply(Eind[,const], MAR=1, all)
    }
    
    if (nE > 0) {
        forbPrimer <- apply(sapply(c(combi$match.F[[x]],combi$match.R[[x]]), FUN=function(y) {!Eind[,y]}), MAR=1, any)
    } else {
        forbPrimer <- FALSE
    }

    
    forbid <- forbType | forbPrimer
    Eind <- Eind[!forbid,]
    jNames <- jNames[!forbid]
    seqs <- seqs[!forbid]    

    ## Determine junctions for inclusion and exclusion reads of each element
    inclJ <- lapply(grep("E", toupper(eName)), FUN=function(y) {
        out <- jNames[Eind[,y]]
        if (any(!grepl(eName[y], out))) {out <- NA}
        unique(out)
    })
    if (nE == 0) {inclJ <- list(jNames)}

    ## Remove duplicated junctions (considering read length)
    junc <- data.frame(name=jNames, seq=seqs, stringsAsFactors=FALSE)
    junc <- junc[!duplicated(junc$name),]
    
    ## If still too short, add the adaptor; trim to required length and test if long enough
    junc$nameFinal <- junc$name
    short <- which(nchar(junc$seq) < readEndUp)
    junc$nameFinal[short] <- paste(junc$name[short], "adaptR", sep="")
    junc$seq[short]  <- paste(junc$seq[short], as.character(reverseComplement(DNAString(opt$adaptR))), sep="")
    junc$seq <- substr(junc$seq, start=opt$trimR1.5 + 1, stop=readEndUp)

    ## Match the junctions used for each element with junction names
    inclJ <- sapply(inclJ, FUN=function(y) {
        out <- which(junc$name %in% y)
        if (length(out) == 0) {
            out <- NA
        } else {
            out <- paste(paste("fw", out, sep=""), collapse=",")
        }
        out
    })

    ## In the case of Alt5 or Alt3, indicate the element that it is attached to
    out.Relative <- rep(NA, length(eName))
    out.Relative[eType == "3"] <- eName[grep("3", eType) + 1]
    out.Relative[eType == "5"] <- eName[grep("5", eType) - 1]
    out.Relative[grep("C", out.Relative)] <- NA
    if (nE > 0) {
        out.Relative <- out.Relative[-c(1, length(eName))]
    } else {
        out.Relative <- NA
    }
    
    out.Event <- eName[-c(1, length(eName))]
    if (length(out.Event) == 0) {out.Event <- NA}

    out.Label <- rep(NA, max(1, nE))
    if (length(grep("E", eName)) > 0) {out.Label[grep("E", eName) - 1] <- paste("A", 1:length(grep("E", eName)), sep="")}
            
    outIncl <- data.frame(Gene        = events$gene[x],
                          Event       = out.Event,
                          Label       = out.Label,
                          Relative    = out.Relative,
                          JunctionsFw = inclJ,
                          stringsAsFactors = FALSE
                          )
            
    ## Return the junctions table and element matching data
    junc$nameFinal <- paste(events$gene[x], "_", "fw", 1:nrow(junc), "_", junc$nameFinal, sep="")
    junc <- junc[,c(3,2)]
    names(junc)[1] <- "name"

    list(junc = junc,
         elem = outIncl
         )    
})


## Generate rev junctions
reverse <- lapply(1:nrow(events), FUN=function(x) {
    elem  <- rev(strsplit(events$structure[x], split="[-:]")[[1]])
    eName <- gsub( "\\[.\\]", "", elem)
    eType <- {
        out <- ifelse(grepl("\\[", elem), "", "A")
        out[grepl("\\[3\\]", elem)] <- "3"
        out[grepl("\\[5\\]", elem)] <- "5"
        out[grepl("\\[i\\]", elem)] <- "I"
        out[grepl("\\[c\\]", elem)] <- "C"
        out[c(1, length(out))]      <- "C"
        out
    }
    nE <- length(combi$E[[x]])
  
    Eind <- matrix(ncol=nE + 2, nrow=2^nE)
    Eind[,c(1, ncol(Eind))] <- TRUE
    if (length(combi$E[[x]]) > 0) {
        for (i in 1:nE) {
            Eind[,i + 1] <- rep(rep(c(FALSE, TRUE), each=2^(i-1)), times=2^nE / (2^i))
        }
    }
    Eind <- Eind[,ncol(Eind):1]
    if (length(dim(Eind)) == 0) {Eind <- matrix(TRUE, nrow=1, ncol=2)}

    lengthsAll <- c(nchar(combi$C2[x]), rev(nchar(combi$E[[x]])), nchar(combi$C1[x]))
    lengths <- lapply(1:nrow(Eind), FUN=function(y) {cumsum(lengthsAll[Eind[y,]])})
    select <- sapply(lengths, FUN=function(y) {  # how many elements needed?
        suppressWarnings(out <- min(which(y >= readEndDn)))
        if (is.infinite(out)) {out <- length(y)}
        out
    })
    
    jNames <- sapply(1:nrow(Eind), FUN=function(y) {  # junction names using only the necessary elements
        paste(eName[which(Eind[y,])[1:select[y]]], collapse="")})
    seqsAll <- c(combi$C2[x], rev(combi$E[[x]]), combi$C1[x])
    seqsAll <- as.character(reverseComplement(DNAStringSet(seqsAll)))
    seqs <- sapply(1:nrow(Eind), FUN=function(y) {   # seqs of the necessary elements
        paste(seqsAll[which(Eind[y,])[1:select[y]]], collapse="")})
      
    ## Remove 'forbidden' combinations due to element type and primer overlap
    forbType <- rep(FALSE, nrow(Eind))
    link5 <- which(eType %in% c("5","I"))
    if (length(link5) > 0) {
        forbType <- forbType | apply(sapply(link5, FUN=function(y) {Eind[,y] & !Eind[,y + 1]}), MAR=1, any)
    }
    link3 <- which(eType %in% c("3","I"))
    if (length(link3) > 0) {
        forbType <- forbType | apply(sapply(link3, FUN=function(y) {Eind[,y] & !Eind[,y - 1]}), MAR=1, any)
    }
    const <- which(eType == "C")
    if (length(const) > 2) {
        forbType <- forbType | !apply(Eind[,const], MAR=1, all)
    }

    if (nE > 0) {
        forbPrimer <- apply(sapply(length(elem) - c(combi$match.F[[x]],combi$match.R[[x]]) + 1, FUN=function(y) {!Eind[,y]}), MAR=1, any)
    } else {
        forbPrimer <- FALSE
    }

    forbid <- forbType | forbPrimer
    Eind <- Eind[!forbid,]
    jNames <- jNames[!forbid]
    seqs <- seqs[!forbid]    

    ## Determine junctions for inclusion and exclusion reads of each element
    inclJ <- lapply(grep("E", toupper(eName)), FUN=function(y) {
        out <- jNames[Eind[,y]]
        if (any(!grepl(eName[y], out))) {out <- NA}
        unique(out)
    })
    if (nE == 0) {inclJ <- list(jNames)}

    ## Remove duplicated junctions (considering read length)
    junc <- data.frame(name=jNames, seq=seqs, stringsAsFactors=FALSE)
    junc <- junc[!duplicated(junc$name),]
    
    ## If still too short, add the adaptor; trim to required length and test if long enough
    junc$nameFinal <- junc$name
    short <- which(nchar(junc$seq) < readEndUp)
    junc$nameFinal[short] <- paste(junc$name[short], "adaptF", sep="")
    junc$seq[short]  <- paste(junc$seq[short], as.character(reverseComplement(DNAString(opt$adaptF))), sep="")
    junc$seq <- substr(junc$seq, start=opt$trimR2.5 + 1, stop=readEndUp)

    ## Match the junctions used for each element with junction names
    inclJ <- sapply(inclJ, FUN=function(y) {
        out <- which(junc$name %in% y)
        if (length(out) == 0) {
            out <- NA
        } else {
            out <- paste(paste("rv", out, sep=""), collapse=",")
        }
        out
    })

    ## In the case of Alt5 or Alt3, indicate the element that it is attached to
    out.Relative <- rep(NA, length(eName))
    out.Relative[eType == "3"] <- eName[grep("3", eType) - 1]
    out.Relative[eType == "5"] <- eName[grep("5", eType) + 1]
    out.Relative[grep("C", out.Relative)] <- NA
    if (nE > 0) {
        out.Relative <- out.Relative[-c(1, length(eName))]
    } else {
        out.Relative <- NA
    }
    
    out.Event <- eName[-c(1, length(eName))]
    if (length(out.Event) == 0) {out.Event <- NA}

    out.Label <- rep(NA, max(1, nE))
    if (length(grep("E", eName)) > 0) {out.Label[grep("E", eName) - 1] <- paste("A", length(grep("E", eName)):1, sep="")}
            
    outIncl <- data.frame(Gene        = events$gene[x],
                          Event       = out.Event,
                          Label       = out.Label,
                          Relative    = out.Relative,
                          JunctionsRv = inclJ,
                          stringsAsFactors = FALSE
                          )[length(out.Event):1,]
            
    ## Return the junctions table and element matching data
    junc$nameFinal <- paste(events$gene[x], "_", "rv", 1:nrow(junc), "_", junc$nameFinal, sep="")
    junc <- junc[,c(3,2)]
    names(junc)[1] <- "name"

    list(junc = junc,
         elem = outIncl
         )    
})


## Generate amplicons - note that there may be more amplicons than junctions because
## reads may not reach deep enough into amplicons to distinguish different ones
amplicons <- lapply(1:nrow(events), FUN=function(x) {
    elem  <- strsplit(events$structure[x], split="[-:]")[[1]]
    eName <- gsub( "\\[.\\]", "", elem)
    eType <- {
        out <- ifelse(grepl("\\[", elem), "", "A")
        out[grepl("\\[3\\]", elem)] <- "3"
        out[grepl("\\[5\\]", elem)] <- "5"
        out[grepl("\\[i\\]", elem)] <- "I"
        out[grepl("\\[c\\]", elem)] <- "C"
        out[c(1, length(out))]      <- "C"
        out
    }
    nE <- length(combi$E[[x]])
  
    Eind <- matrix(ncol=nE + 2, nrow=2^nE)
    Eind[,c(1, ncol(Eind))] <- TRUE  
    if (length(combi$E[[x]]) > 0) {
        for (i in 1:nE) {
            Eind[,i + 1] <- rep(rep(c(FALSE, TRUE), each=2^(i-1)), times=2^nE / (2^i))
        }
    }

    #lengthsAll <- c(nchar(combi$C1[x]), nchar(combi$E[[x]]), nchar(combi$C2[x]))
    #lengths <- lapply(1:nrow(Eind), FUN=function(y) {cumsum(lengthsAll[Eind[y,]])})
    #select <- sapply(lengths, FUN=function(y) {  # how many elements needed?
    #    suppressWarnings(out <- min(which(y >= readEndUp)))
    #    if (is.infinite(out)) {out <- length(y)}
    #    out
    #})
    
    jNames <- sapply(1:nrow(Eind), FUN=function(y) {  # junction names (all elements)
        paste(eName[which(Eind[y,])], collapse="")})
    seqsAll <- c(combi$C1[x], combi$E[[x]],combi$C2[x])
    seqs <- sapply(1:nrow(Eind), FUN=function(y) {  
        paste(seqsAll[which(Eind[y,])], collapse="")})
      
    ## Remove 'forbidden' combinations due to element type and primer overlap
    forbType <- rep(FALSE, nrow(Eind))
    link5 <- which(eType %in% c("5","I"))
    if (length(link5) > 0) {
        forbType <- forbType | apply(sapply(link5, FUN=function(y) {Eind[,y] & !Eind[,y - 1]}), MAR=1, any)
    }
    link3 <- which(eType %in% c("3","I"))
    if (length(link3) > 0) {
        forbType <- forbType | apply(sapply(link3, FUN=function(y) {Eind[,y] & !Eind[,y + 1]}), MAR=1, any)
    }
    const <- which(eType == "C")
    if (length(const) > 2) {
        forbType <- forbType | !apply(Eind[,const], MAR=1, all)
    }
    
    if (nE > 0) {
        forbPrimer <- apply(sapply(c(combi$match.F[[x]],combi$match.R[[x]]), FUN=function(y) {!Eind[,y]}), MAR=1, any)
    } else {
        forbPrimer <- FALSE
    }
    
    forbid <- forbType | forbPrimer
    jNames <- jNames[!forbid]
    seqs <- seqs[!forbid]    
            
    ## Return amplicons and names
    list(name = paste(events$gene[x], jNames, sep="_"),
         elem = seqs
         )    
})


## Generate element table
elemF <- data.frame(Gene        = unlist(lapply(forward, FUN=function(x) {x$elem$Gene})),
                    Event       = unlist(lapply(forward, FUN=function(x) {x$elem$Event})),
                    Label       = unlist(lapply(forward, FUN=function(x) {x$elem$Label})),
                    Relative    = unlist(lapply(forward, FUN=function(x) {x$elem$Relative})),
                    JunctionsFw = unlist(lapply(forward, FUN=function(x) {x$elem$JunctionsFw})),
                    stringsAsFactors = FALSE
                )
elemR <- data.frame(Gene        = unlist(lapply(reverse, FUN=function(x) {x$elem$Gene})),
                    Event       = unlist(lapply(reverse, FUN=function(x) {x$elem$Event})),
                    Label       = unlist(lapply(reverse, FUN=function(x) {x$elem$Label})),
                    Relative    = unlist(lapply(reverse, FUN=function(x) {x$elem$Relative})),
                    JunctionsRv = unlist(lapply(reverse, FUN=function(x) {x$elem$JunctionsRv})),
                    stringsAsFactors = FALSE
                    )
if (!all(paste(elemF$Gene, elemF$Event, elemF$Label, elemF$Relative) ==
         paste(elemR$Gene, elemR$Event, elemR$Label, elemR$Relative))) {
    stop("Something went wrong. Different event elements for fwd and rev.")
}

elem <- cbind(elemF[,1:5], JunctionsRv=elemR[,5])
undet <- is.na(elem$JunctionsFw) & is.na(elem$JunctionsRv)
if (any(undet)) {
    warning("The following elements are undetermined in both directions:\n",
            paste(paste(elem$Gene, elem$Event)[undet], collapse="\n"))
    elem$Gene[undet] <- paste("#", elem$Gene[undet], sep="")
}

write.table(elem[order(elem$Gene, elem$Event),], file=opt$eventOut, row.names=F, col.names=T, quote=F, sep='\t')


## Generate junction and amplicon FASTA files
juncF <- data.frame(name = unlist(lapply(forward, FUN=function(x) {x$junc$name})),
                    seq  = unlist(lapply(forward, FUN=function(x) {x$junc$seq})),
                    stringsAsFactors = FALSE
                )
juncR <- data.frame(name = unlist(lapply(reverse, FUN=function(x) {x$junc$name})),
                    seq  = unlist(lapply(reverse, FUN=function(x) {x$junc$seq})),
                    stringsAsFactors = FALSE
                )


shortF <- opt$readLength - opt$trimR1.5 - opt$trimR1.3 - nchar(juncF$seq)
shortF <- data.frame(name=juncF$name, short=shortF)[order(-shortF),]
if (any(shortF$short > 0)) {
    cat("\nThe following amplicons are shorter than trimmed fwd read:\n")
    cat(paste(paste(shortF$name, shortF$short, sep=": ")[shortF$short > 0], collapse="\n"), "\n\n")
}

shortR <- opt$readLength - opt$trimR2.5 - opt$trimR2.3 - nchar(juncR$seq)
shortR <- data.frame(name=juncR$name, short=shortR)[order(-shortR),]
if (any(shortR$short > 0)) {
    cat("\nThe following amplicons are shorter than trimmed rev read:\n")
    cat(paste(paste(shortR$name, shortR$short, sep=": ")[shortR$short > 0], collapse="\n"), "\n\n")
}

if (any(shortF$short > 0) | any(shortR$short > 0)) {
    cat("Junction libraries were not produced\n")
} else {
    ## Make junction FASTA files
    fa.fw <- character(length=2 * nrow(juncF))
    fa.rv <- character(length=2 * nrow(juncR))
    
    fa.fw[seq(1, length(fa.fw) - 1, 2)] <- paste0(">", juncF$name)
    fa.fw[seq(2, length(fa.fw), 2)]     <- juncF$seq
    
    fa.rv[seq(1, length(fa.rv) - 1, 2)] <- paste0(">", juncR$name)
    fa.rv[seq(2, length(fa.rv), 2)]     <- juncR$seq

    write.table(fa.fw, file=paste(opt$juncOut, "_fwd.fa", sep=""), row.names=F, col.names=F, quote=F, sep='\t')
    write.table(fa.rv, file=paste(opt$juncOut, "_rev.fa", sep=""), row.names=F, col.names=F, quote=F, sep='\t')

    ## Amplicon FASTA file
    fa.amp <- character(2 * length(amplicons))
    fa.amp[seq(1, length(fa.amp) - 1, 2)] <- paste0(">", sapply(amplicons, "[[", 1))
    fa.amp[seq(2, length(fa.amp),     2)] <- sapply(amplicons, "[[", 2)

    write.table(fa.amp, file=paste(opt$juncOut, "_amplicons.fa", sep=""), row.names=F, col.names=F, quote=F, sep='\t')

    ## Produce minimum distance reports
    minDist.fw <- minDistJunc(fa.fw)
    minDist.rv <- minDistJunc(fa.rv)
    if (any(minDist.fw$distAll < 2)) {
        warning(paste(length(which(minDist.fw$distAll < 2)),
                      "fwd junctions have a Hamming distance of < 2 to another fwd junction."))
    } else {
        cat("Minimum Hamming distance of fwd junctions is", min(minDist.fw$distAll), "\n")
    }
    if (any(minDist.rv$distAll < 2)) {
        warning(paste(length(which(minDist.rv$distAll < 2)),
                      "rev junctions have a Hamming distance of < 2 to another rev junction."))
    } else {
        cat("Minimum Hamming distance of rev junctions is", min(minDist.rv$distAll), "\n")
    }
    
    write.table(minDist.fw, file=paste(opt$juncOut, "_fwd_minDist.tab", sep=""), row.names=F, col.names=T, quote=F, sep='\t')
    write.table(minDist.rv, file=paste(opt$juncOut, "_rev_minDist.tab", sep=""), row.names=F, col.names=T, quote=F, sep='\t')

    ## Produce junction BED files to be used by the pipeline
    bed.fw <- data.frame(juncF$name, 0, nchar(juncF$seq))
    bed.rv <- data.frame(juncR$name, 0, nchar(juncR$seq))

    write.table(bed.fw, file=paste(opt$juncOut, "_fwd.bed", sep=""), row.names=F, col.names=F, quote=F, sep='\t')
    write.table(bed.rv, file=paste(opt$juncOut, "_rev.bed", sep=""), row.names=F, col.names=F, quote=F, sep='\t')

}