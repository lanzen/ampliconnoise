#!/usr/bin/perl

use strict;

my $cutoff = $ARGV[0];
my $nOTUs  = 0;
my $minval = 0;

my %OTUVectors = {};

while(my $line = <STDIN>){
    my @tokens = split(/ /,$line);
    my $d = shift(@tokens);

    if($d == $cutoff){
	$nOTUs = shift(@tokens);

	for(my $i = 0; $i < $nOTUs; $i++){
	    my %hashCluster = {};

	    my @Cluster = split(/,/,$tokens[$i]);
	    foreach my $entry(@Cluster){
		$entry =~ /^(\w+)_\d+_(\d+)$/;
		if($1 ne ""){
		    #print "$1 $2\n";
		    $hashCluster{$1} += $2;
		}
	    }

	    foreach my $sample(keys %hashCluster){
		my $ref = $OTUVectors{$sample};
		if($sample =~/^\w+$/){
		    #print "good $sample\n";
		    if($ref ne undef){
			${$OTUVectors{$sample}}[$i] = $hashCluster{$sample}; 
		    }
		    else{
			my @vector = ();
			$vector[$i] = $hashCluster{$sample};
			$OTUVectors{$sample} = \@vector; 
		    }
		}
		else{
		   # print "bad $sample\n";
		}
	    }
	}
    }
}

my @Vectors = ();
my @names = ();
foreach my $sample(sort (keys %OTUVectors)){
    if($OTUVectors{$sample} ne undef){
	my @vector = @{$OTUVectors{$sample}};
		
	for(my $i = 0; $i < $nOTUs; $i++){
	    if($vector[$i] == undef){
		$vector[$i] = 0;
	    }
	}
	
	#print "$sample\n";
	push(@names,$sample);
	push(@Vectors, \@vector);
    }
}

my $nSamples = scalar(@Vectors);

my @print = ();

for(my $j = 0; $j < $nOTUs; $j++){

    $print[$j] = 0;

    for(my $i = 0; $i < $nSamples; $i++){
	if($Vectors[$i][$j] > 6){
	    $print[$j] = 1;
	}
    }

}
printf("OTU,");
my $nstring = join(",",@names);
print "$nstring\n";
for(my $j = 0; $j < $nOTUs; $j++){
    my $i = 0;
    printf("C%d,",$j);
    for($i = 0; $i < $nSamples - 1; $i++){
	printf("%d,", $Vectors[$i][$j]);
    } 
    printf("%d\n",$Vectors[$i][$j]);	
}


