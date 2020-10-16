#-v!/bin/bash
# Kevin Ha and Ulrich Braunschweig, 2013-2020

FILE=$1     # (potentially gzipped) FASTQ file
GENOME=$2   # bowtie index with junction library
OUTDIR=$3   # output directory
CORES=$4
TRIM_5=$5
TRIM_3=$6
BTOPT=$7
BIN=$8

mkdir -vp $OUTDIR

# uncompress file if compressed
INDIR=`dirname $FILE`
FILE=`basename $FILE`
BASE=`basename $FILE .gz`
BASE=`basename $BASE .fq`
BASE=`basename $BASE .fastq`

# Bowtie options
# -t        print the time
# --sam     output in SAM format
# -m        suppress all alignments for a particular read or pair if more than <int>
# reportable alignments exist for it.
# -k        report up to <int> valid alignments per read or pair
# -v        report alignments with at most <int> mismatches
# --un      unmapped reads
# --max     alignments exceeding limit set by -m to <filename>
# --al      reads with at least one alignment report to <filename>

bowtie --sam $BTOPT -p $CORES \
    --trim3 $TRIM_3 --trim5 $TRIM_5 \
    --max $OUTDIR/$BASE.multi.fq \
    --un $OUTDIR/$BASE.nomap.fq \
    $GENOME $INDIR/$FILE | samtools view -b -o $OUTDIR/$BASE.unsorted.bam && \

# Sort BAM file
samtools sort $OUTDIR/$BASE.unsorted.bam -o $OUTDIR/$BASE.bam && \
pigz -vf -9 -p $CORES $OUTDIR/$BASE.*.fq

rm -v $OUTDIR/$BASE.unsorted.bam

$BIN/get_unmapped_stats.sh $OUTDIR/$BASE.bam

echo End: `date`
