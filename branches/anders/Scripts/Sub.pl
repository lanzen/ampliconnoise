#!/usr/bin/perl

$fastaFile = $ARGV[0];

$ucFile = $ARGV[1];

open(FILE, $fastaFile) or die "Can't open $fastaFile\n";
my $count = 0;
my %hashID = {};
while($line = <FILE>){
	chomp($line);

	if($line=~/>(.*)/){
		$hashID{$1} = $count;
		$count++;
	}
}

close(FILE);

open(FILE, $ucFile) or die "Can't open $ucFile\n";

while($line = <FILE>){
        chomp($line);
	
	my @tokens = split(/\t/,$line);
	my $id = $tokens[8];
	my $nid = $hashID{$id};

	$tokens[8] = $nid;
	my $string = join("\t",@tokens);
	print "$string\n";
}

close(FILE);

