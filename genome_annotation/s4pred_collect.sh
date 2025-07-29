#!/bin/bash

s4path=~/data/mammary/processedData/s4pred/
outfile=~/data/mammary/processedData/Homo_sapiens.GRCh38.pep.large_s4pred.fas
 

##rm $outfile -f
for file in $s4path/*
do
    id=`basename $file|sed 's/\\.ss2//'`
    if ! grep -q $id $outfile; then
	##echo adding $id 
	echo ">" $id >> $outfile
	sed 's/^ \+//;s/ \+/\t/g' $file | cut -f 3 | sed -e '1,2d' | tr -d '\n' >> $outfile
	echo "" >> $outfile
    fi
done
