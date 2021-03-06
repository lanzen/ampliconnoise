#!/usr/bin/perl

use POSIX;

my $nTotal = 0;
my $nClean = 0;
my $flowSeq = "TACG";

my $primer    = &translateIUPAC($ARGV[0]);
my $revPrimer = reverse($ARGV[1]);

$revPrimer =~ tr/[A,C,G,T]/[T,G,C,A]/; 
$revPrimer = &translateIUPAC($revPrimer);

my $out       = $ARGV[2];
my $maxFlows = 800;

my $ffile = "${out}.fa";
my $dfile = "${out}.dat";
open(FFILE, ">${ffile}") or die "Can't open ${out}.fa\n";
open(DFILE, ">${dfile}") or die "Can't open ${out}.dat\n";

while(my $line = <STDIN>){
    chomp($line);
 
    if($line =~ />(.*)/){
	my $id = $1;
	#print "$id\n";
	$line = <STDIN>;
	chomp($line);
	#print "$id $line\n";
	my @flows = split(/ /,$line);
	
	#print "@flows\n";
	my $length = shift(@flows);
	my $read;
	my $newLength = 0;
        my $noise  = 0;
	my $signal = 0;
	#print "$length\n";
	#print "$length @flows\n";
	while($newLength*4 < $length){
	    my $signalL = 0;	
	    for($j = 0; $j < 4; $j++){
		my $f = $flows[$j + 4*$newLength];
		if($f > 0.5){
		    $signalL++;
		    if($f < 0.7 || $f > 6.49){
			$noise++;
		    }
		}
		
	    }
    
	    if($signalL == 0){
		$signal++;
	    }

	    $newLength++;
	}
	$newLLength = 4*$newLength;
	#print "$length $newLength $newLLength $noise\n";
	if($noise <=  1 && $signal == 0){
		my $read = flowToSeq($length,@flows);
	#	print "$read $revPrimer\n";	   
		if($read=~/^TCAG.*(${primer}.*${revPrimer}).*$/){
	  	#    	print "Good\n";
	    		printf DFILE "$id $length ";
	    		for($j = 0; $j < $maxFlows; $j++){
				printf DFILE "$flows[$j] ";
	    		}
		
	    		printf DFILE "\n";
		
	    		printf FFILE ">$id\n";
	    		printf FFILE "$1\n";
		
	    		$nClean++;
	    		goto found;
		}

		found : $nTotal++;
    	}
}
}

close(DFILE);
close(FFILE);

open my $in,  '<',  $dfile      or die "Can't read old file: $!";
open my $out, '>', "$dfile.new" or die "Can't write new file: $!";

print $out "$nClean $maxFlows\n";

while( <$in> )
{     
    print $out $_;
}

close $out;

rename("$dfile.new", $dfile);
exit(0);

sub flowToSeq()
{
    my ($length, @flowgram) = @_;
    my $retSeq = "";

    for(my $i = 0; $i < $length; $i++){
	my $signal = floor($flowgram[$i] + 0.5);
	my $base   = substr($flowSeq, $i % 4, 1);

	for(my $j = 0; $j < $signal; $j++){
	    $retSeq .= $base;
	}
    }

    return $retSeq;
}

sub translateIUPAC()
{
  my ($seq) = @_;
  $seq=~s/B/\[CGT\]/g;
  $seq=~s/S/\[GC\]/g;    
  $seq=~s/N/\[ACTG\]/g;
  $seq=~s/Y/\[CT\]/g;
  $seq=~s/R/\[AG\]/g;
  $seq=~s/D/\[AGT\]/g;
  $seq=~s/H/\[ACT\]/g;
  $seq=~s/M/\[AC\]/g;
  $seq=~s/K/\[GT\]/g;
  $seq=~s/V/\[ACG\]/g;
  return $seq;
}
