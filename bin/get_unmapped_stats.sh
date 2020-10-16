#!/bin/bash

for i in $*
do
    B=`basename $i .bam`
    samtools view $i | awk -v NM=$B '{if ($3 == "*") UM = UM + 1} END {OFS="\t"; PCT=100*(NR - UM)/NR; print NM, UM, NR, PCT}' - > `dirname $i`/${B}.stats.txt
done
