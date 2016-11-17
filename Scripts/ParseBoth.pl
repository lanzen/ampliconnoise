#!/usr/bin/perl

my $tag = shift(@ARGV);
my $revTag = shift(@ARGV);

$revTag =~ tr/[A,C,G,T]/[T,G,C,A]/; 
$revTag = &translateIUPAC($revPrimer);

$tag =~ /G*(.*)/;
$tag2 = &translateIUPAC($1);


while($line = <STDIN>){
    chomp($line);

    if($line =~ />(.*)/){
	print "$line\n";
    }
    else{
	$line =~ /$tag2(.*)$revTag.*/;
 	
	$parse = $1;

	print "$parse\n";
    }
}

sub translateIUPAC()
{
  my ($seq) = @_;

  $seq=~s/W/\[AT\]/g;
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


