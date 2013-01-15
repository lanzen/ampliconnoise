#!/bin/bash

export nodes=8
export otu_dist=0.03
export primer=primer.fasta
min_size=50
max_size=2000

#Fixes warning message with uDAPL error message appearing:
export mpiextra="--mca btl tcp,self"

export CLASSPATH=$AMPLICON_NOISE_HOME/lib/ampliconflow.jar:$AMPLICON_NOISE_HOME/lib/core-1.8.1.jar
export PYRO_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/LookUp_Titanium.dat
export SEQ_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/Tran.dat 

spyro=60
#PyroNoise cluster size

cpyro=0.01
#PyroNoise cluster init


if test $# -le 0; then
   echo "Usage: RunPreSplitXL.sh filename.sff [primersequence]"
   exit -1
fi

if [ -n "$2" ]; then
	primer=$2
elif [ ! -f primer.fasta ]; then
	echo "Can't find file primer.fasta containing the primer sequence!"
	exit
fi

stub=${1//.sff}

# first generate sff text file if necessary                                                                   
if [ ! -f ${stub}.sff.txt ]; then
    echo "Generating .sff.txt file"
    sffinfo $1 >${stub}.sff.txt
fi  

echo "Parsing sff.txt file"
java ampliconflow.sff.FlowsFlex $stub.sff.txt $primer

export cropFL=`tail -1 $stub.stat.txt`

#Unique sequences - large datasets only
echo "Getting unique sequences"
FastaUnique -in ${stub}.filtered.fasta > ${stub}_U.fa


#sequence distances - large datasets only
echo "Calculating sequence distances"
usearch -cluster_fast ${stub}_U.fa -id 0.70 -centroids ${stub}_U_c.fasta -uc ${stub}_U.uc
Sub.pl ${stub}_U.fa ${stub}_U.uc > ${stub}_U.ucn
if [ ! -d ${stub}_split ]; then
    mkdir ${stub}_split
    cp ${stub}.dat ${stub}.filtered.map ${stub}_U.ucn ${stub}_split
fi
cd ${stub}_split
	
SplitClusterClust -din ${stub}.dat -min ${stub}.filtered.map -uin ${stub}_U.ucn -m 100 > ${stub}_split.stats


#mpirun $mpiextra -np $nodes NDist -i -in ${stub}_U.fa > ${stub}_U_I.ndist
xs=$?
if [[ $xs != 0 ]]; then
    echo "Presplitting exited with status $xs"
    exit $xs
fi


#Run PyroDist and PyroNoise on each cluster separetely
echo "Calculating .fdist files using PyroDist"
for c in C*
do
        if [ -d $c ] ; then
            mpirun $mpiextra -np $nodes PyroDist -in ${c}/${c}.dat -out ${c}/${c} > ${c}/${c}.fout
	    xs=$?
	    if [[ $xs != 0 ]]; then
		echo "PyroDist exited with status $xs"
		exit $xs
	    fi
        fi
done

echo "Clustering .fdist files using FCluster"

for c in C*
do
        if [ -d $c ] ; then
	    FCluster -in ${c}/${c}.fdist -out ${c}/${c}_X > ${c}/${c}.fout
	    xs=$?
	    if [[ $xs != 0 ]]; then
		echo "FCluster exited with status $xs"
		exit $xs
	    fi
	    rm ${c}/${c}.fdist
        fi
done

echo "Running PyroNoiseM Clustering"
for dir in C*
do
        if [ -d $dir ] ; then
            mpirun $mpiextra -np $nodes PyroNoiseM -din ${dir}/${dir}.dat -out ${dir}/${dir}_s60_c01 -lin ${dir}/${dir}_X.list -s 60.0 -c 0.01 > ${dir}/${dir}_s60_c01.pout
	    xs=$?
	    if [[ $xs != 0 ]]; then
		echo "PyroNoiseM exited with status $xs"
		exit $xs
	    fi
        fi
done

cat C*/C*_s60_c01_cd.fa > ${stub}_s60_c01_cd.fa
cat C*/C*_s60_c01.mapping > ${stub}_s60_c01.mapping

echo "Cropping barcodes, primes and low quality end (at 400 bp)"
Truncate.pl 400 < ${stub}_s60_c01_cd.fa > ${stub}_s60_c01_T400.fa
cropF.py  ${stub}_s60_c01_T400.fa $cropFL > ${stub}_s60_c01_T400_P_BC.fa

echo "Running SeqDist"
mpirun $mpiextra -np $nodes SeqDist -in ${stub}_s60_c01_T400_P_BC.fa > ${stub}_s60_c01_T400_P_BC.seqdist
xs=$?
if [[ $xs != 0 ]]; then
    echo "SeqDist exited with status $xs"
    exit $xs
fi


echo "Clustering SeqDist output"
FCluster -in ${stub}_s60_c01_T400_P_BC.seqdist -out ${stub}_s60_c01_T400_P_BC_S > ${stub}_s60_c01_T400_P_BC.fcout
xs=$?
if [[ $xs != 0 ]]; then
    echo "FCluster exited with status $xs"
    exit $xs
fi

echo "Running SeqNoise"
mpirun $mpiextra -np $nodes SeqNoise -in ${stub}_s60_c01_T400_P_BC.fa -din ${stub}_s60_c01_T400_P_BC.seqdist -lin ${stub}_s60_c01_T400_P_BC_S.list -out ${stub}_s60_c01_T400_P_BC_s30_c08 -s 30.0 -c 0.08 -min ${stub}_s60_c01.mapping > ${stub}_s60_c01_T400_P_BC_s30_c08.snout
xs=$?
if [[ $xs != 0 ]]; then
    echo "SeqNoise exited with status $xs"
    exit $xs
fi
rm *.seqdist

ln -s ${stub}_s60_c01_T400_P_BC_s30_c08_cd.fa ${stub}_F.fa

echo "Running Perseus"
Perseus -sin ${stub}_F.fa > ${stub}_F.per
xs=$?
if [[ $xs != 0 ]]; then
    echo "Persus exited with status $xs"
    exit $xs
fi

Class.pl ${stub}_F.per -7.5 0.5 > ${stub}_F.class
FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa
rm Temp.*

echo "Clustering OTUs"
mpirun $mpiextra -np $nodes NDist -i -in ${stub}_F_Good.fa > ${stub}_F_Good.ndist
xs=$?
if [[ $xs != 0 ]]; then
    echo "NDist exited with status $xs"
    exit $xs
fi

FCluster  -i -in ${stub}_F_Good.ndist -out ${stub}_F_Good > ${stub}_F_Good.fdist
xs=$?
if [[ $xs != 0 ]]; then
    echo "FCluster exited with status $xs"
    exit $xs
fi

echo "Writing otu representatives, calculating Rarefaction and Chao1 estimates"

java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -repseq ${stub}_F_Good.fa > ${stub}_OTUs_${otu_dist}.fasta

us=`grep -ce ">" ${stub}_F.fa`
cs=`grep -ce ">" ${stub}_F_Chi.fa`
otus=`grep -ce ">" ${stub}_OTUs_${otu_dist}.fasta`

cd ..
ln -s ${stub}_split/${stub}_F_Good.fa .
ln -s ${stub}_split/${stub}_F_Good.list .


if [ ! -f AN_stats.txt ]; then
    echo -e 'Sample\tTotal reads\tPre-filtered reads\tUnique sequences\tChimeric sequences\tRemaining unique sequences\tRemaining reads\tOTUs\tShannon index\tSimpsons index (1-D)' > AN_stats.txt
fi

tr=`grep -ce ">" ${stub}.raw.fasta`
pf=`grep -ce ">" ${stub}.filtered.fasta`
rus=`grep -ce ">" ${stub}_F_Good.fa`
rr=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -s -weigh -totalreads`
shannon=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -weigh -s -shannon`
simpson=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -weigh -s -simpson`

echo -e "${stub}\t${tr}\t${pf}\t${us}\t${cs}\t${rus}\t${rr}\t${otus}\t${shannon}\t${simpson}" >> AN_stats.txt
