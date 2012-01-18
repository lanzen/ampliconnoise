export bc=keys.csv
export nodes=4
export otu_dist=0.03

#Fixes warning message with uDAPL error message appearing:
export mpiextra="--mca btl tcp,self" 

export PYRO_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/LookUp_E123.dat
export SEQ_LOOKUP_FILE=$AMPLICON_NOISE_HOME/Data/Tran.dat 

if test $# -le 0; then
   echo "Usage: RunFLX.sh filename.sff [primersequence]"
   exit -1
fi

if [ -n "$2" ]; then
	primer=$2
elif [ ! -f primer.fasta ]; then
	echo "Can't find file primer.fasta containing the primer sequence!"
	exit
else
	while read line; do    
    	if [ "${line:0:1}" != ">" ]; then
    		primer=$line
    		break
    	fi
done <primer.fasta   
fi
echo "Primer sequence: $primer"

stub=${1//.sff}
stub=${stub//.txt}

# first generate sff text file if necessary                                                                   
if [ ! -f ${stub}.sff.txt ]; then
    echo "Generating .sff.txt file"
    sffinfo $1 >${stub}.sff.txt
fi

pLength=`expr ${#primer}`

echo "Parsing sff.txt file"
if [ -f $bc ]; then
    echo "Using barcodes file $bc for splitting"
    SplitKeys.pl $primer $bc < ${stub}.sff.txt > matching.fasta 2>nonmatching.fasta
    firstKey=`head -1 $bc`
    firstKey=${firstKey//*,}
    bcLength=${#firstKey}
else
    echo "No barcodes file found. Using entire dataset"
    FlowsFA360.pl $primer $stub< ${stub}.sff.txt
    touch ${stub}.raw
    bcLength=0
fi

cropFL=`expr $bcLength + $pLength`

if [ ! -f AN_stats.txt ]; then
    echo -e 'Sample\tTotal reads\tPre-filtered reads\tUnique sequences\tChimeric sequences\tRemaining unique sequences\tRemaining reads\tOTUs\tShannon index\tSimpsons index (1-D)' > AN_stats.txt
fi

for file in *.raw; do
    stub=${file//.raw}
    if [ ! -f ${stub}.dat ]; then
	Clean360.pl $primer $stub < $file
    fi

    echo "Running PyroDist for ${stub}"
    mpirun $mpiextra -np $nodes PyroDist -in ${stub}.dat -out ${stub} > ${stub}.fout

    echo "Clustering PyroDist output for ${stub}"
    FCluster -in ${stub}.fdist -out ${stub}_X > ${stub}.fout

    rm ${stub}.fdist
    rm ${stub}_X.otu ${stub}_X.tree

    echo "Running PyronoiseM for ${stub}"
    mpirun $mpiextra -np $nodes PyroNoiseM -din ${stub}.dat -out ${stub}_s60_c01 -lin ${stub}_X.list -s 60.0 -c 0.01 > ${stub}_s60_c01.pout

    echo "Cropping barcodes, primes and low quality end (at 220 bp)"
    Truncate.pl 220 < ${stub}_s60_c01_cd.fa > ${stub}_s60_c01_T220.fa
    cropF.py  ${stub}_s60_c01_T220.fa $cropFL > ${stub}_s60_c01_T220_P_BC.fa

    echo "Running SeqDist for ${stub}"
    mpirun $mpiextra -np $nodes SeqDist -in ${stub}_s60_c01_T220_P_BC.fa > ${stub}_s60_c01_T220_P_BC.seqdist

    echo "Clustering SeqDist output for ${stub}"
    FCluster -in ${stub}_s60_c01_T220_P_BC.seqdist -out ${stub}_s60_c01_T220_P_BC_S > ${stub}_s60_c01_T220_P_BC.fcout

    echo "Running SeqNoise for ${stub}"
    mpirun $mpiextra -np $nodes SeqNoise -in ${stub}_s60_c01_T220_P_BC.fa -din ${stub}_s60_c01_T220_P_BC.seqdist -lin ${stub}_s60_c01_T220_P_BC_S.list -out ${stub}_s60_c01_T220_P_BC_s30_c08 -s 30.0 -c 0.08 -min ${stub}_s60_c01.mapping > ${stub}_s60_c01_T220_P_BC_s30_c08.snout

    ln -s ${stub}_s60_c01_T220_P_BC_s30_c08_cd.fa ${stub}_F.fa

    echo "Running Perseus for ${stub}"
    Perseus -sin ${stub}_F.fa > ${stub}_F.per
    Class.pl ${stub}_F.per -7.5 0.5 > ${stub}_F.class
    FilterGoodClass.pl ${stub}_F.fa ${stub}_F.class 0.5 1>${stub}_F_Chi.fa 2>${stub}_F_Good.fa

    echo "Clustering OTUs for ${stub}"
    mpirun $mpiextra -np $nodes NDist -i -in ${stub}_F_Good.fa > ${stub}_F_Good.ndist

    FCluster -i -in ${stub}_F_Good.ndist -out ${stub}_F_Good > ${stub}_F_Good.fdist

    echo "Writing otu representatives "

    java amliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -repseq ${stub}_F_Good.fa > ${stub}_OTUs_${otu_dist}.fasta

    tr=`grep -ce ">" ${stub}.sff.txt`
    pf=`head -1 ${stub}.dat`
    pf=${pf//" "*}
    us=`grep -ce ">" ${stub}_F.fa`
    cs=`grep -ce ">" ${stub}_F_Chi.fa`
    rus=`grep -ce ">" ${stub}_F_Good.fa`
    rr=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -s -weigh -totalreads`
    otus=`grep -ce ">" ${stub}_OTUs_${otu_dist}.fasta`
    shannon=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -weigh -s -shannon`
    simpson=`java ampliconflow.otu.OTUUtils -in ${stub}_F_Good.list -dist $otu_dist -weigh -s -simpson`

    echo -e "${stub}\t${tr}\t${pf}\t${us}\t${cs}\t${rus}\t${rr}\t${otus}\t${shannon}\t${simpson}" >> AN_stats.txt


done

