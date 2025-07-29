#!/bin/sh


bioawk=/home/${USER}/programs/bioawk/
iupred=/home/${USER}/programs/iupred3/
outdir=/home/${USER}/data/mammary/processedData/iupred3/


idfile=/home/${USER}/data/mammary/processedData/Homo_sapiens.GRCh38.pep.large_ids.txt
fasta=/home/${USER}/data/mammary/processedData/Homo_sapiens.GRCh38.pep.large.fa

iutype=long
smooth=medium
min_len=20

tmpfa=tmp.fa
for id in $(cat $idfile); do
    echo $id
    outfile=${outdir}/${id}_iupred3.tsv
    if [ ! -f $outfile ]; then
	${bioawk}/bioawk -v var="${id}" -c fastx \
		 '$name ~ var {print ">"$name"\n"$seq}' $fasta > ${tmpfa}
	len=`${bioawk}/bioawk -c fastx '{ print length($seq) }' ${tmpfa}`
	if [[ $len -ge $min_len ]]; then
	    python3 ${iupred}/iupred3.py -a -s $smooth ${tmpfa}  $iutype \
		| awk 'NR == 1, /# POS/ { next } { print }' > $outfile
	else
	    echo -e "\t" too short: ${len}.
	fi
    else
	echo -e "\t" exists.
    fi
done

## MANUAL STEP: do this loop in a different terminal in parallel
## to above

## do the same but comming from bottom up,
## using a different fasta file name
pmtfa=pmt.fa
for id in $(tac $idfile); do
    echo $id
    outfile=${outdir}/${id}_iupred3.tsv
    if [ ! -f $outfile ]; then
	${bioawk}/bioawk -v var="${id}" -c fastx \
		 '$name ~ var {print ">"$name"\n"$seq}' $fasta > ${pmtfa}
	len=`${bioawk}/bioawk -c fastx '{ print length($seq) }' ${pmtfa}`
	if [[ $len -ge $min_len ]]; then
	    python3 ${iupred}/iupred3.py -a -s $smooth ${pmtfa}  $iutype \
		| awk 'NR == 1, /# POS/ { next } { print }' > $outfile
	else
	    echo -e "\t" too short: ${len}.
	fi
    else
	echo -e "\t" exists.
    fi
done

