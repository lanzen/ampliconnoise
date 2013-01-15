#!/usr/bin/perl

use strict;

my $mapFile = $ARGV[0];

open(MAPFILE, $mapFile) or die "Can't open $mapFile\n";

my @map = ();

while(my $line = <MAPFILE>){
    chomp($line);
    my @tokens = split(/:/,$line);
    push(@map, $tokens[2]);
}

close(MAPFILE);


my $nOTUs  = 0;

while(my $line = <STDIN>){
    chomp($line);
    my @tokens = split(/ /,$line);
    my $d = shift(@tokens);

    $nOTUs = shift(@tokens);
    
    print "$d $nOTUs ";

    for(my $i = 0; $i < $nOTUs; $i++){
	my @tokens2 = split(/,/,$tokens[$i]);
	my $size = scalar(@tokens2);

	for(my $j = 0; $j < $size - 1; $j++){
	    print "$map[$tokens2[$j]],";
	}
	print "$map[$tokens2[$size - 1]] ";
    }
    print "\n";
}


sub min
{
    my ($x, $y) = @_;

    if($x < $y){
	return $x;
    }
    else{
	return $y;
    }
}
