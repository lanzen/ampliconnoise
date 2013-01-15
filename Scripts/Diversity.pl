#!/usr/bin/perl

$line = <STDIN>;
chomp($line);

@samples = split(/,/,$line);
my @vector = ();
shift(@samples);

$nSamples = scalar(@samples);
$j = 0;
while($line = <STDIN>){
	chomp($line);

	my @tokens = split(/,/,$line);

	shift(@tokens);

	for($i = 0; $i < $nSamples; $i++){
		$vector[$i][$j] = $tokens[$i];
	}
	$j++; 
}
$nOTUs = $j;

my @total = ();
my @S = ();
my @chao = ();
my @shannon = ();
my @P = ();

for($i = 0; $i < $nSamples; $i++){
	my $n1 = 0;
	my $n2 = 0;

	for($j = 0; $j < $nOTUs; $j++){
		$total[$i]+=$vector[$i][$j];
		if($vector[$i][$j] > 0){
			$S[$i]++;
		}
		if($vector[$i][$j] == 1){
			$n1++;
		}
		if($vector[$i][$j] == 2){
			$n2++;
		}
	}

	for($j = 0; $j < $nOTUs; $j++){
		if($vector[$i][$j] > 0){
			$p = $vector[$i][$j] / $total[$i];
			$shannon[$i] += -$p*log($p);
		}
	}
	$P[$i] = $shannon[$i]/log($S[$i]); 
	if($n1 > 0 && $n2 > 0){
		$chao[$i] = $S[$i] + 0.5*(($n1*$n1)/$n2);
	}
}
printf("%s,%s\n","Diversity",join(',',@samples));
printf("%s,%s\n","N",join(',',@total));
printf("%s,%s\n","S",join(',',@S));

print "Chao,";

for($i = 0; $i < $nSamples; $i++){
	printf("%.2f,",$chao[$i]);
}
printf("\n");
print "Shannon,";

for($i = 0; $i < $nSamples; $i++){
	printf("%.2f,",$shannon[$i]);
}
printf("\n");
print "Pielou,";

for($i = 0; $i < $nSamples; $i++){
	printf("%.2f,",$P[$i]);
}
printf("\n");
