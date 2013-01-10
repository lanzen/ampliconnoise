#!/bin/bash

inputdir=$1
burn=$2

printf "Sample, Log-normal, Invese Gaussian, Sichel\n";
for file in $inputdir/*_LN_0.sample
do
	stub=${file%_LN_0.sample};
	stub=${stub#$inputdir/}
	printf "$stub,"
	MeanSeries.pl 4 $burn ${inputdir}/${stub}_LN_*.sample
	MeanSeries.pl 4 $burn ${inputdir}/${stub}_IG_*.sample
	MeanSeries.pl 5 $burn ${inputdir}/${stub}_SI_*.sample
	printf "\n";
done
