#!/bin/bash

defaultSize=500

newsize=${2:-$defaultSize} #second argument primer as a Perl regular expression

for datfile in *.dat
do
	stub=${datfile%.dat};
	stub=${stub#${sampledir}};
	echo $datfile $stub $newsize

	#generate flowgram and fasta files
	SubsampleDat.pl $datfile $newsize > ${stub}_S${newsize}.dat
done

