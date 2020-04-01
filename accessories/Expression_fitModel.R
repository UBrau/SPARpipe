#!/usr/bin/env Rscript
### Differential expression analysis of SPAR-seq screen
### Read table containing reads per event and sample, generate a design matrix
### based on a treatment table, and estimate the dispersion in an edgeR GLM model.
### Save model for the second step (extracting a specific contrast).
###
### U. Braunschweig, 10/2017


libMissing <- !require(optparse, quietly=T) && stop("Failed to load R package 'optparse'")

### Parse input
args <- commandArgs(TRUE)
option.list <- list(
    make_option(c("-r", "--readsTab"),       action="store", type="character", metavar="FILE",
                help="Path of table containing reads per event in each sample from the 3_combine.pl
                module of SPARpipe"),
    make_option(c("-t", "--treatTab"),       action="store", type="character", metavar="FILE",
                help="Path of treatment table used in previous SPARpipe modules"),
    make_option(c("-d", "--designTab"),      action="store", type="character", metavar="FILE",
                help="Path of design factor table where lines are the same as in the treatment table
                and columns correspond to factors that can be added to the model.
                Required columns are 'ID' and 'Replicate' like in the other tables; at least one
                column should be present representing treatment, with negative controls representing
                the baseline having the same name"),
    make_option(c("-e", "--designElem"),     action="store", type="character", metavar="STRING",
                help="Which factors to model in the design matrix. Contrasts must correspond to
                column names in the design factor table. Additionally, the reference state for each
                factor must be specified, in the form: Factor1:State1,Factor2:State2,..."),
    make_option(c("-c", "--cores"),          default=1,
                help="CPU cores to run with [%default]")
)
parser <- OptionParser(option_list=option.list,
                       usage="usage: %prog [options] DIR",
                       description="Model dispersion and fit GLM from counts table with edgeR")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) != 1)
    stop(print_help(parser))


#### Function definitions ####

simplifyCounts <- function(x, const=c("Replicate","ID")) {
    out <- as.matrix(x[,!(names(x) %in% const)])
    rownames(out) <- paste(x$Replicate, x$ID, sep=".")
    genes <- sub("\\.[eE][0-9]+$", "", colnames(out))

    duplos <- unique(genes[duplicated(genes)])
    if (length(duplos) > 0) {
        duploInd <- sapply(duplos, FUN=function(y) {
            ind <- as.integer(names(sort(table(apply(out[,genes ==  y], MAR=1, which.max)), decreasing=T)))
            which(genes == y)[ind]
        })
    }

    ind <- sort(c(which(genes %in% setdiff(genes, duplos)), duploInd))
    out <- out[,ind]
    colnames(out) <- genes[ind]
    out
}
    

#### Main part ####

### Load packages
libMissing <- !require(edgeR, quietly=T) && stop("Failed to load R package 'edgeR'")


### Check input
inDir <- sub("\\/*$","",opt$args[1])

if (is.null(opt$options$readsTab))               {stop("--readsTab must be specified")}
if (is.null(opt$options$treatTab))               {stop("--treatTab must be specified")}
if (is.null(opt$options$designTab))              {stop("--designTab must be specified")}
if (is.null(opt$options$designElem))             {stop("--designElem must be specified")}

if (!dir.exists(inDir))                          {stop("Input directory not found")}
if (!dir.exists(file.path(inDir, "raw")))        {stop("Subdirectory raw/ not found")}
if (!dir.exists(file.path(inDir, "expression"))) {dir.create(file.path(inDir, "expression"))}
if (!file.exists(opt$options$readsTab))          {stop("Read count table not found")}
if (!file.exists(opt$options$treatTab))          {stop("Treatment table not found")}
if (!file.exists(opt$options$designTab))         {stop("Design factor table not found")}
if (!file.exists(file.path(inDir, "raw", "ReadsPerEvent.tab"))) {
    stop("Read count table not found at ", paste(inDir, "/raw/ReadsPerEvent.tab", sep=""))
}

desf <- strsplit(strsplit(opt$options$designElem, split=",")[[1]], split=":")
desf <- data.frame(factor = sapply(desf, "[", 1),
                   base   = sapply(desf, "[", 2),
                   stringsAsFactors = F)

## Load tables and check content
reads  <- read.delim(opt$options$readsTab, check.names=F)
treat  <- read.delim(opt$options$treatTab)
destab <- read.delim(opt$options$designTab)

if (!all(paste(reads$Replicate, reads$ID) == paste(treat$Replicate, treat$ID))) {
    stop("Read and treatment tables could not be matched based on columns 'Replicate' and 'ID'")
}
if (!all(paste(reads$Replicate, reads$ID) == paste(destab$Replicate, destab$ID))) {
    stop("Read and contrast tables could not be matched based on columns 'Replicate' and 'ID'")
}

if (!all(desf$factor %in% names(destab))) {
    stop("Columns not found in design factor table: ",
         paste(desf$factor[!(desf$factor %in% names(design))], collapse = ", "))
}
for (i in 1:nrow(desf)) {
    if (!(desf$base[i] %in% destab[,desf$factor[i]])) {
        stop("No '", desf$base[i], "' treatments found for design factor ", desf$factor[i])
    }
}


## Prepare counts matrix and design matrix
counts <- simplifyCounts(reads)
write.table(data.frame(Replicate = reads$Replicate,
                       ID        = reads$ID,
                       counts,
                       check.names=F),
            row.names=F, col.names=T, quote=F, sep='\t',
            file=file.path(inDir, "expression", "ReadsPerGene.tab")
            )

designfac <- lapply(desf$factor, FUN=function(x) {
    relevel(as.factor(destab[,x]), desf$base[desf$factor == x])
})
names(designfac) <- desf$factor

desform <- formula(paste("~ ", paste(desf$factor, collapse=" + "), sep=""))
design <- model.matrix(desform, data=model.frame(desform, data=designfac))


## Estimate dispersion
dge <- DGEList(counts = t(counts))
dge <- calcNormFactors(dge)

system.time(dge <- estimateDisp(dge, design=design))
save(dge, file=file.path(inDir, "expression", "edgeR.estDispersion.Robj"))

system.time(fit <- glmQLFit(dge, design=design, robust=T))
save(dge, fit, file=file.path(inDir, "expression", "edgeR.modelFit.Robj"))
file.remove(file.path(inDir, "expression", "edgeR.estDispersion.Robj"))

cat("R archive containing fitted model saved\n")