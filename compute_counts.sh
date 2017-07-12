#!/bin/bash

for i in $*
do
    BASE=`basename $i _Fwd.sorted.bam`
    DIR=`dirname $i`
    DIR=`dirname $DIR`
    UP=/home/blencowe/blencowe31/ulrich/proj/SPARseq.HongAndy.neural.170202/seq/junctions/Junctions108_11-101_fwd.bed
    DOWN=/home/blencowe/blencowe31/ulrich/proj/SPARseq.HongAndy.neural.170202/seq/junctions/Junctions108_11-101_rev.bed
    echo $BASE

    bedtools coverage -counts -b $i \
        -a $UP \
        | cut -f 1,4 - \
        | awk -F "\t" '{OFS="\t"; split($1,n,"_"); print n[1], n[2], n[3], $2}' - \
        | sed 's/up/up_/' - \
        | sed 's/dn/dn_/' - \
        > ${DIR}/coverage/${BASE}.cov.tab

    bedtools coverage -counts -b `dirname $1`/${BASE}_Rev.sorted.bam \
        -a $DOWN \
        | cut -f 1,4 - \
        | awk -F "\t" '{OFS="\t"; split($1,n,"_"); print n[1], n[2], n[3], $2}' - \
        | sed 's/up/up_/' - \
        | sed 's/dn/dn_/' - \
        >> ${DIR}/coverage/${BASE}.cov.tab

    #./3_compute_psi_rpm_values.R coverage/${BASE}.cov.tab
done
