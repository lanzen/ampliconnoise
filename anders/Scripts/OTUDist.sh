#!/bin/bash
stub=${1//.fasta}
stub=${stub//.fa}

export CLASSPATH=$AMPLICON_NOISE_HOME/lib/ampliconflow.jar:$AMPLICON_NOISE_HOME/lib/core-1.8.1.jar

nodes=4
export mpiextra="--mca btl tcp,self"

sed -i -e 's/\.00//' $1

FastaUnique -in $1 > $stub.U.fa

mpirun $mpiextra -np $nodes NDist -i -in $stub.U.fa > ${stub}.U.ndist
FCluster -i -in ${stub}.U.ndist -out $stub.U > $stub.U.fdist

java ampliconflow.otu.OTUUtils -in $stub.U.list -addfrom $stub.map -weigh -distable > $stub.otudist.txt
java ampliconflow.otu.OTUUtils -in $stub.U.list -addfrom $stub.map -repseq $1 > $stub.otus.fasta
