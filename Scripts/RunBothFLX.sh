#!/bin/bash

#barcode file
bc=keys.csv
nodes=32
snodes=1
min_size=50
max_size=40000
#Fixes warning message with uDAPL error message appearing:
mpiextra="--mca btl tcp,self" 

export PYRO_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/LookUp_E123.dat
export SEQ_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/Tran.dat

#hardcoded parameters for AmpliconNoise

length=400
#truncation length

spyro=60
#PyroNoise cluster size

cpyro=0.01
#PyroNoise cluster init

sseq=30
#SeqNoise cluster size

cseq=0.08
#SeqNoise cluster init

alpha=-7.5
#Perseus logit intercept

beta=0.5
#Perseus logit gradient

#file locations

primerfile=primer.fasta
reversefile=reverse.fasta

#read in primer sequence
if [ ! -f $primerfile ]; then
        echo "Can't find file $primerfile containing the primer sequence!"
        exit
else
        while read line; do
        if [ "${line:0:1}" != ">" ]; then
                primer=$line
                break
        fi
        done < $primerfile
fi
echo "Primer sequence: $primer"

#read in primer sequence
if [ ! -f $reversefile ]; then
        echo "Can't find file $reversefile containing the reverse primer sequence!"
        exit
else
        while read line; do
        if [ "${line:0:1}" != ">" ]; then
                reverse=$line
                break
        fi
        done < $reversefile
fi
echo "Primer sequence: $reverse"



split()
{
	echo "split:"
# first generate sff text file if necessary
        echo $stub 	                                                           
	if [ ! -f ${stub}.sff.txt ]; then
    		echo "Generating .sff.txt file"
    		sffinfo ${stub}.sff >${stub}.sff.txt
	fi

	echo "Parsing sff.txt file"
	if [ -f $bc ]; then
		if [ -f ${stub}.sff.txt ]; then
    			echo "Using barcodes file $bc for splitting"
    			echo "$primer $bc"
			SplitKeysBoth.pl $primer $reverse $bc < ${stub}.sff.txt > matching.fasta 2>nonmatching.fasta
		fi
	else
    		echo "No barcode file found aborting..."
	fi
}

if [ -f AN_stats.txt ]; then
	rm AN_stats.txt
fi

if [ -f All_Good.fa ]; then
	rm All_Good.fa
fi

echo -e 'Sample\tTotal reads\tPre-filtered reads\tUnique sequences\tChimeric sequences\tRemaining unique sequences\tRemaining clean reads' > AN_stats.txt

filter()
{
	echo "filter:"
	while IFS=, read stub barcode
	do 
    		file=${stub}.raw
    		if [ -f ${stub}.raw ]; then	
			CleanBoth.pl $primer $reverse $stub < $file
    		fi
	done < $bc
}

case $1 in
	pyronoise)
        ;;
        seqnoise)
        ;;
        perseus)
        ;;
        otus)
        ;;
	split)
        if [ $# -eq 2 ]; then
                stub=${2//.sff}
                stub=${stub//.txt}
                split
		exit
        fi
	;;
	filter)
	filter
	exit
	;;
        all)
	if [ $# -eq 2 ]; then
		stub=${2//.sff}
        	stub=${stub//.txt}
		split
        	filter
	fi
	;;
	*)
        	echo "Usage: RunTitanium.sh all sfffile"
        	echo "RunTitanium.sh [pyronoise|seqnoise|perseus|otus|filter|all]"
        	exit
        ;;

esac

#read in primer sequence
if [ ! -f $primerfile ]; then
        echo "Can't find file $primerfile containing the primer sequence!"
        exit
else
        while read line; do
        if [ "${line:0:1}" != ">" ]; then
                primer=$line                break
        fi
        done < $primerfile
fi
echo "Primer sequence: $primer"

count=0

while IFS=, read stub barcode revBarcode
do
                barcodes[$count]=$barcode
                stubs[$count]=$stub

                cleanSizes[$count]=$(head -n 1 ${stub}.dat | cut -d" " -f1)
                let count=count+1
done < $bc



pyronoise()
{
	if [ ! -f ${pstub}_cd.fa ]; then
        	echo "Running PyroDist for ${stub}"

        	mpirun $mpiextra -np $nodes PyroDist -in ${stub}.dat -out ${stub} > ${stub}.fout
        	xs=$?
       	 	if [[ $xs != 0 ]]; then
                	echo "PyroDist exited with status $xs"
                	exit $xs
        	fi


        	echo "Clustering PyroDist output for ${stub}"

        	FCluster -in ${stub}.fdist -out ${stub}_X > ${stub}.fout
        	xs=$?
        	if [[ $xs != 0 ]]; then
                	echo "FCluster exited with status $xs"
                	exit $xs
        	fi

        	rm ${stub}.fdist ${stub}_X.otu ${stub}_X.tree

        	echo "Running PyronoiseM for ${stub}"
        	mpirun $mpiextra -np $nodes PyroNoiseM -din ${stub}.dat -out ${pstub} -lin ${stub}_X.list -s $spyro -c $cpyro > ${pstub}.pout
        	xs=$?
        	if [[ $xs != 0 ]]; then
                	echo "PyroNoiseM parsing exited with status $xs"
                	exit $xs
        	fi

        	echo "Cropping barcodes, primes and low quality end (at 400 bp)"
	fi
}

pyronoisesplit()
{
	echo "Splitting ${stub}"
	
	#get unique sequences
	#if [ ! -f ${stub}_U.fa ]; then
    		echo "Getting unique sequences"
    		FastaUnique -in ${stub}.fa > ${stub}_U.fa
	#fi

	#use NDist to get sequence distances 
	if [ ! -f ${stub}_U.uc ]; then
    		echo "Clustering with uclust"
		usearch -cluster_fast ${stub}_U.fa -id 0.75 -centroids ${stub}_U_c.fasta -uc ${stub}_U.uc
		Sub.pl ${stub}_U.fa ${stub}_U.uc > ${stub}_U.ucn
	fi
	if [ ! -d ${stub}_split ]; then
		mkdir ${stub}_split
		cp ${stub}.dat ${stub}.map ${stub}_U.ucn ${stub}_split
	fi
	cd ${stub}_split

	
	SplitClusterClust -din ${stub}.dat -min ${stub}.map -uin ${stub}_U.ucn -m 100 > ${stub}_split.stats
	

	echo "Calculating .fdist files"
	for c in C*
	do
        	if [ ! -f ${c}/${c}.fdist ] ; then
                	mpirun -np $nodes PyroDist -in ${c}/${c}.dat -out ${c}/${c} > ${c}/${c}.fout
        	fi
	done

	echo "Clustering .fdist files"

	for c in C*
	do
        	if [ ! -f ${c}/${c}_X.list ] ; then
                	FCluster -in ${c}/${c}.fdist -out ${c}/${c}_X > ${c}/${c}.fout
			rm ${c}/${c}.fdist
        	fi
	done

	echo "Running PyroNoise"
	for dir in C*
	do
        	if [ ! -f ${dir}/${dir}_s${spyro}_cd.fa ] ; then
                	mpirun -np $nodes PyroNoiseM -din ${dir}/${dir}.dat -out ${dir}/${dir}_s${spyro} -lin ${dir}/${dir}_X.list -s $spyro -c $cpyro > ${dir}/${dir}_${spyro}.pout
        	fi
	done

	cat C*/C*_s${spyro}_cd.fa > All_s${spyro}_cd.fa
	cat C*/C*_s${spyro}.mapping > All_s${spyro}.mapping
        
	echo "Cropping barcodes, primes and low quality end (at 400 bp)"
	
	sed "s/>.*\(_[0-9]\+_[0-9]\+\)/>${pstub}\1/" All_s${spyro}_cd.fa > ../${pstub}_cd.fa
	
	cp All_s${spyro}.mapping ../${pstub}.mapping
	cd ..
}


seqnoise()
{
	if [ ! -f ${sstub}_cd.fa ]; then
		scount=`grep -ce ">" ${pstub}_cd.fa`
		if [ $scount -lt 100 ]; then 
			tnodes=$snodes
		else
			tnodes=$nodes
		fi

		echo "Running SeqDist for ${stub}"
        	mpirun $mpiextra -np $tnodes SeqDist -in ${pstub}_B.fa > ${pstub}_B.seqdist
        	xs=$?
        	if [[ $xs != 0 ]]; then
        		echo "SeqDist exited with status $xs"
                	exit $xs
        	fi

        	echo "Clustering SeqDist output for ${stub}"
        	FCluster -in ${pstub}_B.seqdist -out ${pstub}_B > ${pstub}_B.fcout
        	xs=$?
        	if [[ $xs != 0 ]]; then
        		echo "FCluster exited with status $xs"
                	exit $xs
        	fi

        	echo "Running SeqNoise for ${stub}"
       		 mpirun $mpiextra -np $tnodes SeqNoise -in ${pstub}_B.fa -din ${pstub}_B.seqdist -lin ${pstub}_B.list -out ${sstub} -s $sseq -c $cseq -min ${pstub}.mapping > ${sstub}.snout
        	xs=$?
        	if [[ $xs != 0 ]]; then
        		echo "SeqNoise exited with status $xs"
                	exit $xs
        	fi
	fi
}

perseus()
{
        echo "Running PerseusD for ${stub}"
        #del=s${spyro}_T${length}_s${sseq}_
        #sed "s/$del//g" ${sstub}_cd.fa > ${stub}_F.fa
	#if [ ! -f ${stub}_F.per ]; then
	#	mkdir ${stub}_F_Per

	#	cp ${stub}_F.fa ${stub}_F_Per
		
	#	cd ${stub}_F_Per
        	
	#	Perseus2 -sin ${stub}_F.fa > ${stub}_F.per&

	#	cd ..
	#fi

	if [ ! -f ${stub}_F_Good.fa ]; then
        	Class.pl ${stub}_F_Per/${stub}_F.per ${alpha} ${beta} > ${stub}_F.class
        	xs=$?
        	if [[ $xs != 0 ]]; then
                	echo "Persus exited with status $xs"
                	exit $xs
		fi
        	./FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa
	fi
}


perseusD()
{
	echo "Running PerseusD for ${stub}"
        del=s${spyro}_T${length}_s${sseq}_
        sed "s/$del//g" ${sstub}_cd.fa > ${stub}_F.fa

        PerseusD -sin ${stub}_F.fa > ${stub}_F.class
        xs=$?
        if [[ $xs != 0 ]]; then
        	echo "Persus exited with status $xs"
                exit $xs
        fi

        FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa
}

otus()
{
	echo "Now construct OTUs across whole data set"
	FastaUnique -in All_Good.fa > All_Good_U.fa
	mpirun $mpiextra -np $nodes NDist -in All_Good_U.fa > All_Good_U.ndist
	FCluster -a -r 0.001 -in All_Good_U.ndist -out All_Good_U
	Map.pl All_Good.map < All_Good_U.list > All_Good.list
	cut -d" " -f1,2 All_Good.list > All_Good.plot

	CSV.pl 0.01 < All_Good.list > All_Good_C01.csv
	CSV.pl 0.03 < All_Good.list > All_Good_C03.csv
	CSV.pl 0.05 < All_Good.list > All_Good_C05.csv
	CSV.pl 0.10 < All_Good.list > All_Good_C10.csv

	Typical.pl 0.01 All_Good.fa All_Good.list > All_Good_C01.fa
	Typical.pl 0.03 All_Good.fa All_Good.list > All_Good_C03.fa
	Typical.pl 0.05 All_Good.fa All_Good.list > All_Good_C05.fa
	Typical.pl 0.10 All_Good.fa All_Good.list > All_Good_C10.fa

	Diversity.pl < All_Good_C01.csv > All_Good_C01.d
	Diversity.pl < All_Good_C03.csv > All_Good_C03.d
	Diversity.pl < All_Good_C05.csv > All_Good_C05.d
	Diversity.pl < All_Good_C10.csv > All_Good_C10.d
}

i=0
while [ $i -lt $count ]
do
	echo "Sample = ${stubs[$i]}"
	stub=${stubs[$i]}
	barcode=${barcodes[$i]}
	size=${cleanSizes[$i]}
	pstub=${stub}_s${spyro}
	sstub=${pstub}_B_s${sseq}
	echo "Number of reads = ${size}"
	if [ "$size" -gt "$min_size" ] ; then
		case $1 in 
		pyronoise)
			if [ ! -f ${pstub}_cd.fa ]; then 
				if [ "$size" -lt "$max_size" ] ; then
					pyronoise  
				else
					pyronoisesplit
				fi
			fi	
			;;
		seqnoise)
			echo "${barcode}${primer} ${reverse} < ${pstub}_cd.fa"
			ParseBoth.pl ${barcode}${primer} ${reverse} < ${pstub}_cd.fa > ${pstub}_B.fa
			seqnoise				
			;;
		perseus)
			if [ ! -f ${stub}_F_Good.fa ]; then
				perseus	
			fi	
			;;
		otus)
			;;
		split)
			;;
		filter)
			;;
		all)

			if [ "$size" -lt "$max_size" ] ; then
                        	pyronoise
                        else
                                pyronoisesplit
                        fi

			seqnoise
			perseus
			;;
		esac

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

case $1 in
	pyronoise)
        ;;
        seqnoise)
        ;;
        perseus)
        ;;
	filter)
	;;
        otus)
        otus               
	;;
        all)
	otus
        ;;
esac
