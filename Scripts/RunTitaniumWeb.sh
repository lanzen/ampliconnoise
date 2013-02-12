#!/bin/bash

#Variables needed for mpirun and PyroNoise 
#(change these settings when installing webapp on new system)
export MPI_HOME=/usr/lib64/openmpi/1.4-gcc

#Max no. of unique sequences per sample after PyroNoise
#(Prevents very long running times caused by chimeras and other artefacts)
SEQ_LIMIT=5000

#number of cores to run on
nodes=12
snodes=2

#barode file
bc=keys.csv
lastline=$(tail -n 1 keys.csv; echo x); lastline=${lastline%x}
if [ "${lastline: -1}" != $'\n' ]; then
    echo >> keys.csv
fi

min_size=50
max_size=2000

#Fixes warning message with uDAPL error message appearing:
mpiextra="--mca btl tcp,self" 

export PATH=$MPI_HOME/bin:$PATH

export PYRO_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/LookUp_Titanium.dat
export SEQ_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/Tran.dat

#hardcoded parameters for AmpliconNoise

length=400
#truncation length

spyro=60
#PyroNoise cluster size

cpyro=0.01
#PyroNoise cluster init

sseq=25
#SeqNoise cluster size

cseq=0.08
#SeqNoise cluster init

alpha=-7.5
#Perseus logit intercept

beta=0.5
#Perseus logit gradient

minflows=360

maxflows=720

#file locations

primerfile=primer.fasta
lastline=$(tail -n 1 primer.fasta; echo x); lastline=${lastline%x}
if [ "${lastline: -1}" != $'\n' ]; then
    echo >> primer.fasta
fi

otudist=$2

stub=${1//.sff}
stub=${stub//.txt}

if [ -f AN_Progress.txt ]; then
	rm AN_Progress.txt
fi

if [ -f AN_stats.txt ]; then
	rm AN_stats.txt
fi
echo -e 'Sample\tTotal\tPre-filtered\tUnique\tChimeric\tCleanSeq\tCleanReads' > AN_stats.txt

if [ -f All_Good.fa ]; then
	rm All_Good.fa
fi

#read in primer sequence

while read line; do
    if [ "${line:0:1}" != ">" ]; then
        primer=$line
        break
    fi
done < $primerfile

xs=$?
if [[ $xs != 0 ]]; then
    echo "Error: Could not read primer" >> AN_Progress.txt
    exit $xs
fi

echo "Primer sequence: $primer"  >> AN_Progress.txt

echo "Generating .sff.txt file"  >> AN_Progress.txt
sffinfo $1 >${stub}.sff.txt
xs=$?
if [[ $xs != 0 ]]; then
    echo "Error: SFF parsing exited with status $xs" >> AN_Progress.txt
    exit $xs
fi

echo "Splitting barcodes..."  >> AN_Progress.txt

SplitKeys.pl $primer $bc < ${stub}.sff.txt > matching.fasta 2>nonmatching.fasta

xs=$?
if [[ $xs != 0 ]]; then
    echo "Error: Could not split barcodes. SplitKeys exited with status $xs" >> AN_Progress.txt
    exit $xs
fi

echo "Filtering..." >> AN_Progress.txt

while IFS=, read stub barcode
do 
    file=${stub}.raw
    if [ -f ${stub}.raw ]; then	
	CleanMinMax.pl $primer $stub $minflows $maxflows < $file
    fi
done < $bc

xs=$?
if [[ $xs != 0 ]]; then
    echo "Error: Could not split barcodes. CleanMinMax.pl exited with status $xs" >> AN_Progress.txt
    exit $xs
fi

count=0

while IFS=, read stub barcode
do
    barcodes[$count]=$barcode
    stubs[$count]=$stub
    
    cleanSizes[$count]=$(head -n 1 ${stub}.dat | cut -d" " -f1)
    let count=count+1
done < $bc

i=0
while [ $i -lt $count ]
do
    stub=${stubs[$i]}
    echo "Sample = $stub" >> AN_Progress.txt
    barcode=${barcodes[$i]}
    size=${cleanSizes[$i]}
    pstub=${stub}_s${spyro}
    sstub=${pstub}_T${length}_s${sseq}
    echo "Number of reads = ${size}" >> AN_Progress.txt

    if [ $size -lt $min_size ] ; then
	echo "WARNING: Insufficient reads remain after filtering. The dataset will not be filtered" >> AN_Progress.txt
    elif [ $size" -lt $max_size ] ; then

	## ___ RUNNING NORMAL PYRONOISE AND SEQNOISE ___
	echo "Running PyroDist for ${stub}" >> AN_Progress.txt
	
	mpirun $mpiextra -np $nodes PyroDist -in ${stub}.dat -out ${stub} > ${stub}.fout
	xs=$?
	if [[ $xs != 0 ]]; then
	    echo "PyroDist exited with status $xs"  >> AN_Progress.txt
	    exit $xs
	fi
	echo "Clustering PyroDist output for ${stub}"  >> AN_Progress.txt
	
	FCluster -in ${stub}.fdist -out ${stub}_X > ${stub}.fout
	xs=$?
	if [[ $xs != 0 ]]; then
	    echo "FCluster exited with status $xs"  >> AN_Progress.txt
	    exit $xs
	fi
	
	rm ${stub}.fdist ${stub}_X.otu ${stub}_X.tree
	
	echo "Running PyroNoise for ${stub}"  >> AN_Progress.txt
	mpirun $mpiextra -np $nodes PyroNoiseM -din ${stub}.dat -out ${pstub} -lin ${stub}_X.list -s $spyro -c $cpyro > ${pstub}.pout
        xs=$?
        if [[ $xs != 0 ]]; then
            echo "PyroNoiseM parsing exited with status $xs"  >> AN_Progress.txt
            exit $xs
        fi
	
        echo "Cropping barcodes, primes and low quality end (at 400 bp)"  >> AN_Progress.txt
        Parse.pl ${barcode}${primer} $length < ${pstub}_cd.fa > ${pstub}_T${length}.fa
	
	if [ ! -f ${sstub}_cd.fa ]; then
	    scount=`grep -ce ">" ${pstub}_cd.fa`
	    echo "Counted $scount unique sequences after PyroNoise" >> AN_Progress.txt

	    if [ $scount -lt 100 ]; then 
		tnodes=$snodes
	    else
		tnodes=$nodes
	    fi
	    
	    if [ $scount -gt $SEQ_LIMIT ]; then
		echo "ABORTING RUN: Sample $stub contains too many sequences after PyroNoise step."  >> AN_Progress.txt
		echo "Due to time restrictions your job is therefore cancelled. Please contact the service group for assistance (services@bioinfo.no)"  >> AN_Progress.txt
		exit 134
	    fi
	else
	    echo "Error: $sstub_cd.fa does not exist" >> AN_Progress.txt
	    exit -1
	fi
	
	echo "Running SeqDist for ${stub}"  >> AN_Progress.txt 
        mpirun $mpiextra -np $tnodes SeqDist -in ${pstub}_T${length}.fa > ${pstub}_T${length}.seqdist
        xs=$?
        if [[ $xs != 0 ]]; then
            echo "Error: SeqDist exited with status $xs" >> AN_Progress.txt
            exit $xs
        fi

        echo "Clustering SeqDist output for ${stub}"  >> AN_Progress.txt
        FCluster -in ${pstub}_T${length}.seqdist -out ${pstub}_T${length} > ${pstub}_T${length}.fcout
        xs=$?
        if [[ $xs != 0 ]]; then
            echo "Error: FCluster exited with status $xs"  >> AN_Progress.txt
            exit $xs
        fi
	echo "Running SeqNoise for ${stub}"  >> AN_Progress.txt
       	mpirun $mpiextra -np $tnodes SeqNoise -in ${pstub}_T${length}.fa -din ${pstub}_T${length}.seqdist -lin ${pstub}_T${length}.list -out ${sstub} -s $sseq -c $cseq -min ${pstub}.mapping > ${sstub}.snout
        xs=$?
        if [[ $xs != 0 ]]; then
            echo "Error: SeqNoise exited with status $xs"  >> AN_Progress.txt
            exit $xs
        fi

	echo "Running PerseusD for ${stub}" >> AN_Progress.txt
	del=s${spyro}_T${length}_s${sseq}_
	sed "s/$del//g" ${sstub}_cd.fa > ${stub}_F.fa
	
	PerseusD -sin ${stub}_F.fa > ${stub}_F.class
	xs=$?
	if [[ $xs != 0 ]]; then
            echo "Error: Persus exited with status $xs" >> AN_Progress.txt
            exit $xs
	fi
	
	FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa
	
    else
	
	    ## ___ RUNNING PRESPLITTING FOLLOWED BY PYRONOISE AND SEQNOISE ___
	echo "Splitting ${stub}"  >> AN_Progress.txt
	
            #get unique sequences
    	echo "Getting unique sequences"  >> AN_Progress.txt
    	FastaUnique -in ${stub}.fa > ${stub}_U.fa
	
             #use usearch to get sequence distances 
	echo "Clustering with usearch"  >> AN_Progress.txt
	usearch -cluster_fast ${stub}_U.fa -id 0.70 -centroids ${stub}_U_c.fasta -uc ${stub}_U.uc > /dev/null
	Sub.pl ${stub}_U.fa ${stub}_U.uc > ${stub}_U.ucn
	
	if [ ! -d ${stub}_split ]; then
	    mkdir ${stub}_split
	    cp ${stub}.dat ${stub}.map ${stub}_U.ucn ${stub}_split
	fi
	cd ${stub}_split
	
	SplitClusterClust -din ${stub}.dat -min ${stub}.map -uin ${stub}_U.ucn -m 100 > ${stub}_split.stats

	echo "Running PyroDist (clustered mode)"  >> AN_Progress.txt
	for c in C*
	do
            mpirun $mpiextra -np $nodes PyroDist -in ${c}/${c}.dat -out ${c}/${c} > ${c}/${c}.fout
	done
	
	echo "Clustering PyroDist output (clustered mode)"  >> AN_Progress.txt
	
	for c in C*
	do
            FCluster -in ${c}/${c}.fdist -out ${c}/${c}_X > ${c}/${c}.fout
	    rm ${c}/${c}.fdist
	done
	
	echo "Running PyroNoise (clustered mode)"  >> AN_Progress.txt
	for dir in C*
	do
            mpirun $mpiextra -np $nodes PyroNoiseM -din ${dir}/${dir}.dat -out ${dir}/${dir}_s${spyro} -lin ${dir}/${dir}_X.list -s $spyro -c $cpyro > ${dir}/${dir}_${spyro}.pout
	done
	
	echo "Cropping barcodes, primes and low quality end (at 400 bp; clustered mode)"  >> AN_Progress.txt
	
	for dir in C*
	do
	    Parse.pl ${barcode}${primer} $length < ${dir}/${dir}_s${spyro}_cd.fa > ${dir}/${dir}_s${spyro}_T${length}.fa
	done
	
	echo "Concatenating noise-cleaned sequenecs" >> AN_Progress.txt
	cat C*/C*_s${spyro}_cd.fa > All_s${spyro}_cd.fa
	cat C*/C*_s${spyro}.mapping > All_s${spyro}.mapping

	controlcount=`grep -ce ">" All_s${spyro}_cd.fa`
	echo "Counted $controlcount unique sequences after PyroNoise" >> AN_Progress.txt
	
	if [ $controlcount -gt $ SEQ_LIMIT ]; then
	    echo "ABORTING RUN! Sample $stub contains too many sequences after PyroNoise step."  >> AN_Progress.txt
	    echo "Due to time restrictions your job is therefore cancelled. Please contact the service group for assistance (services@bioinfo.no)"  >> AN_Progress.txt
	    exit 134
	fi

	sed "s/>.*\(_[0-9]\+_[0-9]\+\)/>${pstub}\1/" All_s${spyro}_cd.fa > ../${pstub}_cd.fa
	
	cp All_s${spyro}.mapping ../${pstub}.mapping
	cd ..
	Parse.pl ${barcode}${primer} $length < ${pstub}_cd.fa > ${pstub}_T${length}.fa
	
        cd ${stub}_split
	echo "Running SeqDist and SeqNoise (Clustered mode)"  >> AN_Progress.txt
	for dir in C*
        do
	    cd $dir
	    
	    Parse.pl ${barcode}${primer} $length < ${dir}_s${spyro}_cd.fa > ${dir}_s${spyro}_T${length}.fa
            mpirun $mpiextra -np $nodes SeqDist -in ${dir}_s${spyro}_T${length}.fa > ${dir}_s${spyro}_T${length}.seqdist
            xs=$?
            if [[ $xs != 0 ]]; then
                echo "Error: SeqDist exited with status $xs"  >> AN_Progress.txt
                exit $xs
            fi
	    
            FCluster -in ${dir}_s${spyro}_T${length}.seqdist -out ${dir}_s${spyro}_T${length} > ${dir}_s${spyro}_T${length}.fcout
	    
	    xs=$?
            if [[ $xs != 0 ]]; then
                echo "Error: FCluster exited with status $xs"  >> AN_Progress.txt
                exit $xs
            fi
	    
            mpirun $mpiextra -np $nodes SeqNoise -in ${dir}_s${spyro}_T${length}.fa -din ${dir}_s${spyro}_T${length}.seqdist -lin ${dir}_s${spyro}_T${length}.list -out ${dir}_s${spyro}_T${length}_s${sseq} -s $sseq -c $cseq -min ${dir}_s${spyro}.mapping > ${dir}_s${spyro}.snout
            xs=$?
            if [[ $xs != 0 ]]; then
                echo "Error: SeqNoise exited with status $xs"  >> AN_Progress.txt
                exit $xs
            fi
	    cd ..
	done
	
	cat C*/C*_s${spyro}_T${length}_s${sseq}_cd.fa > ../${sstub}_A.fa
        cat C*/C*_s${spyro}_T${length}_s${sseq}_cd.mapping > ../${sstub}_A.mapping
	
	cd ..
	
	echo "Running SeqDist for All"  >> AN_Progress.txt
	
        scount=`grep -ce ">" ${sstub}_A.fa`
	
        mpirun $mpiextra -np $nodes SeqDist -in ${sstub}_A.fa > ${sstub}_A.seqdist
        
	xs=$?
        
	if [[ $xs != 0 ]]; then
            echo "Error: SeqDist exited with status $xs"  >> AN_Progress.txt
            exit $xs
        fi
	
        echo "Clustering SeqDist output for All"  >> AN_Progress.txt
        
	FCluster -in ${sstub}_A.seqdist -out ${sstub}_A > ${sstub}_A.fcout
	
        xs=$?
        if [[ $xs != 0 ]]; then
            echo "Error: FCluster exited with status $xs"  >> AN_Progress.txt
            exit $xs
        fi
	
        echo "Running SeqNoise for All"  >> AN_Progress.txt
        
	mpirun $mpiextra -np $nodes SeqNoise -in ${sstub}_A.fa -din ${sstub}_A.seqdist -lin ${sstub}_A.list -out ${sstub} -s $sseq -c $cseq -min ${sstub}_A.mapping > ${sstub}_A.snout
        
	xs=$?
        if [[ $xs != 0 ]]; then
            echo "Error: SeqNoise exited with status $xs"  >> AN_Progress.txt
            exit $xs
        fi

	echo "Running PerseusD for ${stub}" >> AN_Progress.txt
	del=s${spyro}_T${length}_s${sseq}_
	sed "s/$del//g" ${sstub}_cd.fa > ${stub}_F.fa
	
	PerseusD -sin ${stub}_F.fa > ${stub}_F.class
	xs=$?
	if [[ $xs != 0 ]]; then
            echo "Error: Persus exited with status $xs" >> AN_Progress.txt
            exit $xs
	fi
	
	FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa
	
    fi
    
    if [ -f ${stub}.raw ] ; then
        tr=`grep -ce ">" ${stub}.raw`
    else
        tr=na
    fi
    
    if [ -f ${stub}.dat ]; then
        pf=`head -1 ${stub}.dat`
        pf=${pf//" "*}
    else
        pf=na
    fi
    
    if [ -f ${stub}_F.fa ]; then
        us=`grep -ce ">" ${stub}_F.fa`
    else
        us=na
    fi
    
    if [ -f ${stub}_F_Chi.fa ]; then
        cs=`grep -ce ">" ${stub}_F_Chi.fa`
    else
        cs=na
    fi
    
    if [ -f ${stub}_F_Good.fa ]; then
        rus=`grep -ce ">" ${stub}_F_Good.fa`
        ccread=`CountFasta.pl < ${stub}_F_Good.fa`
        cat ${stub}_F_Good.fa >> All_Good.fa
    else
        rus=na
        ccread=na
    fi
    
    echo -e "${stub}\t${tr}\t${pf}\t${us}\t${cs}\t${rus}\t$ccread" >> AN_stats.txt
    let i=i+1
done

cat ${stub}_F_Good.fa >> All_Good.fa

echo "Constructing OTUs across whole data set" >> AN_Progress.txt
FastaUnique -in All_Good.fa > All_Good_U.fa
mpirun $mpiextra -np $nodes NDist -in All_Good_U.fa > All_Good_U.ndist
FCluster -a -r 0.001 -in All_Good_U.ndist -out All_Good_U
Map.pl All_Good.map < All_Good_U.list > All_Good.list
cut -d" " -f1,2 All_Good.list > All_Good.plot

CSV.pl $otudist < All_Good.list > OTU_table.csv
Typical.pl $otudist All_Good.fa All_Good.list > OTUs.fa
Diversity.pl < OTU_table.csv > Diversity.csv

