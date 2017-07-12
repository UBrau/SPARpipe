#-v!/bin/bash
# Kevin Ha and Ulrich Braunschweig, 2013-2017

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
unpigz -vf -p $CORES $FILE
INDIR=`dirname $FILE`
FILE=`basename $FILE .gz`
BASE=`basename $FILE .fq`
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
cat $INDIR/$FILE | bowtie -t --sam $BTOPT -p $CORES \
    --trim3 $TRIM_3 --trim5 $TRIM_5 \
    --max $OUTDIR/$BASE.multi.fq \
    --un $OUTDIR/$BASE.nomap.fq \
    $GENOME - | samtools view -b -o $OUTDIR/$BASE.bam && \

# recompress file
pigz -vf -9 -p $CORES $INDIR/$FILE

# Sort BAM file
samtools sort $OUTDIR/$BASE.bam -o $OUTDIR/$BASE.sorted.bam && \

pigz -vf -9 -p $CORES $OUTDIR/$BASE.*.fq

rm -v $OUTDIR/$BASE.bam

$BIN/get_unmapped_stats.sh $OUTDIR/$BASE.sorted.bam

echo End: `date`
