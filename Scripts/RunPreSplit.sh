export primer=primer.fasta
export nodes=4
export otu_dist=0.03

#Fixes warning message with uDAPL error message appearing:
export mpiextra="--mca btl tcp,self"


export CLASSPATH=$AMPLICON_NOISE_HOME/lib/ampliconflow.jar:$AMPLICON_NOISE_HOME/lib/core-1.8.1.jar
export PYRO_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/LookUp_Titanium.dat
export SEQ_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data//Data/Tran.dat 

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

echo "Running PyroDist"
mpirun $mpiextra -np $nodes PyroDist -in ${stub}.dat -out ${stub} > ${stub}.fout

echo "Clustering PyroDist output"
FCluster -in ${stub}.fdist -out ${stub}_X > ${stub}.fout

echo "Running Pyronoise"
mpirun $mpiextra -np $nodes PyroNoiseM -din ${stub}.dat -out ${stub}_s60_c01 -lin ${stub}_X.list -s 60.0 -c 0.01 > ${stub}_s60_c01.pout

echo "Cropping barcodes, primes and low quality end (at 400 bp)"
Truncate.pl 400 < ${stub}_s60_c01_cd.fa > ${stub}_s60_c01_T400.fa
cropF.py  ${stub}_s60_c01_T400.fa $cropFL > ${stub}_s60_c01_T400_P_BC.fa

echo "Running SeqDist"
mpirun $mpiextra -np $nodes SeqDistM -in ${stub}_s60_c01_T400_P_BC.fa > ${stub}_s60_c01_T400_P_BC.seqdist


echo "Clustering SeqDist output"
FClusterM -in ${stub}_s60_c01_T400_P_BC.seqdist -out ${stub}_s60_c01_T400_P_BC_S > ${stub}_s60_c01_T400_P_BC.fcout

echo "Running SeqNoise"
mpirun $mpiextra -np $nodes SeqNoiseM -in ${stub}_s60_c01_T400_P_BC.fa -din ${stub}_s60_c01_T400_P_BC.seqdist -lin ${stub}_s60_c01_T400_P_BC_S.list -out ${stub}_s60_c01_T400_P_BC_s30_c08 -s 30.0 -c 0.08 -min ${stub}_s60_c01.mapping > ${stub}_s60_c01_T400_P_BC_s30_c08.snout

ln -s ${stub}_s60_c01_T400_P_BC_s30_c08_cd.fa ${stub}_F.fa

echo "Running Perseus"
Perseus -sin ${stub}_F.fa > ${stub}_F.per
Class.pl ${stub}_F.per -7.5 0.5 > ${stub}_F.class
FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa

echo "Clustering OTUs"
mpirun $mpiextra -np $nodes NDist -i -in ${stub}_F_Good.fa > ${stub}_F_Good.ndist

FCluster -i -in ${stub}_F_Good.ndist -out ${stub}_F_Good > ${stub}_F_Good.fdist

echo "Writing otu representatives"

java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -repseq ${stub}_F_Good.fa > ${stub}_OTUs_${otu_dist}.fasta

#RankAbundance.pl $otu_dist < ${stub}_F_Good.list > ${stub}_${otu_dist}.adist
#Chao.pl < ${stub}_${otu_dist}.adist > ${stub}_${otu_dist}.chao
#ERarefaction -in ${stub}_${otu_dist}.adist -out ${stub}_${otu_dist} > ${stub}_${otu_dist}.rare

