#!/bin/bash

inputdir=$1
burn=$2
effort=$3

printf "Sample, Log-normal, Invese Gaussian, Sichel\n";
for file in $inputdir/*_LN_0.sample
do
	stub=${file%_LN_0.sample};
	stub=${stub#$inputdir/}
	printf "$stub,"

	cat ${inputdir}/${stub}_LN_*.sample > ${inputdir}/${stub}_LN.catsample
	cat ${inputdir}/${stub}_IG_*.sample > ${inputdir}/${stub}_IG.catsample  
	cat ${inputdir}/${stub}_SI_*.sample > ${inputdir}/${stub}_SI.catsample

	LNRarefaction -in ${inputdir}/${stub}_LN.catsample -a ${stub}.adist -b $burn -s 100 -c $effort 
	printf ",";
	IGRarefaction -in ${inputdir}/${stub}_IG.catsample -a ${stub}.adist -b $burn -s 100 -c $effort 
	
	printf ",";
	SIRarefaction -in ${inputdir}/${stub}_SI.catsample -a ${stub}.adist -b $burn -s 100 -c $effort 
	
	printf "\n";
done
