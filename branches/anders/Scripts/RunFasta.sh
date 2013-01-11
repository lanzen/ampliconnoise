export nodes=2
export otu_dist=0.03
export primer=primer.fasta

#Fixes warning message with uDAPL error message appearing:
export mpiextra="--mca btl tcp,self"

export CLASSPATH=$AMPLICON_NOISE_HOME/lib/ampliconflow.jar:$AMPLICON_NOISE_HOME/lib/core-1.8.1.jar

if test $# -le 0; then
   echo "Usage: RunFasta.sh filename.fasta"
   exit -1
fi

if [ -n "$2" ]; then
	primer=$2
elif [ ! -f primer.fasta ]; then
	echo "Can't find file primer.fasta containing the primer sequence!"
	exit
fi

stub=${1//.fasta}
stub=${stub//.fa}

#FastaUnique -in $1 > ${stub}_U.fa

echo "Cropping barcodes, primes and low quality end (at 400 bp)"
Truncate.pl 400 < ${stub}.fa > ${stub}_T400.fa
#cropF.py  ${stub}_T400.fa $cropFL > ${stub}_T400.fa

echo "Running SeqDist"
mpirun $mpiextra -np $nodes SeqDist -in ${stub}_T400.fa > ${stub}_T400.seqdist
xs=$?
if [[ $xs != 0 ]]; then
    echo "SeqDist exited with status $xs"
    exit $xs
fi

echo "Clustering SeqDist output"
FCluster -in ${stub}_T400.seqdist -out ${stub}_T400_S > ${stub}_T400.fcout
xs=$?
if [[ $xs != 0 ]]; then
    echo "FCluster exited with status $xs"
    exit $xs
fi

echo "Running SeqNoise"
mpirun $mpiextra -np $nodes SeqNoise -in ${stub}_T400.fa -din ${stub}_T400.seqdist -lin ${stub}_T400_S.list -out ${stub}_T400_s30_c08 -s 30.0 -c 0.08 > ${stub}_T400_s30_c08.snout
xs=$?
if [[ $xs != 0 ]]; then
    echo "SeqNoise exited with status $xs"
    exit $xs
fi

ln -sf ${stub}_T400_s30_c08_cd.fa ${stub}_F.fa

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

FCluster -i -in ${stub}_F_Good.ndist -out ${stub}_F_Good > ${stub}_F_Good.fdist
xs=$?
if [[ $xs != 0 ]]; then
    echo "FCluster exited with status $xs"
    exit $xs
fi

echo "Writing otu representatives and statistics"

java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -repseq ${stub}_F_Good.fa > ${stub}_OTUs_${otu_dist}.fasta

if [ ! -f AN_stats.txt ]; then
    echo -e 'Sample\tTotal reads\tUnique sequences\tChimeric sequences\tRemaining unique sequences\tRemaining reads\tOTUs\tShannon index\tSimpsons index (1-D)' > AN_stats.txt
fi

tr=`grep -ce ">" $1`
us=`grep -ce ">" ${stub}_F.fa`
cs=`grep -ce ">" ${stub}_F_Chi.fa`
rus=`grep -ce ">" ${stub}_F_Good.fa`
rr=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -s -weigh -totalreads`
otus=`grep -ce ">" ${stub}_OTUs_${otu_dist}.fasta`
shannon=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -weigh -s -shannon`
simpson=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -weigh -s -simpson`

echo -e "${stub}\t${tr}\t${pf}\t${us}\t${cs}\t${rus}\t${rr}\t${otus}\t${shannon}\t${simpson}" >> AN_stats.txt
