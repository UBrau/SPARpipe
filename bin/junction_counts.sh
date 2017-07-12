#!/bin/bash

FILE=$1
JUNC=$2

BASE=`basename $FILE _Fwd.sorted.bam`
#BASE=`basename $BASE _R1.sorted.bam`
DIR=`dirname $FILE`
DIR=`dirname $DIR`

echo "Calculating junction counts for $BASE"

bedtools coverage -counts -a ${JUNC}_fwd.bed -b $FILE \
    | cut -f 1,4 - \
    | awk -F "\t" '{OFS="\t"; split($1,n,"_"); print n[1], n[2], n[3], $2}' - \
    | sed 's/up/up_/' - \
    | sed 's/dn/dn_/' - \
    > ${DIR}/counts/${BASE}.counts.tab

#bedtools coverage -counts -a ${JUNC}_rev.bed -b `dirname $FILE`/${BASE}_R2.sorted.bam \
bedtools coverage -counts -a ${JUNC}_rev.bed -b `dirname $FILE`/${BASE}_Rev.sorted.bam \
    | cut -f 1,4 - \
    | awk -F "\t" '{OFS="\t"; split($1,n,"_"); print n[1], n[2], n[3], $2}' - \
    | sed 's/up/up_/' - \
    | sed 's/dn/dn_/' - \
    >> ${DIR}/counts/${BASE}.counts.tab


