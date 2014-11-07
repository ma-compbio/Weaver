#!/usr/bin/perl
#
#
#open(I,"<",shift@ARGV);
$CUTOFF = shift@ARGV || 8;
#print $CUTOFF;
while(<STDIN>){
	chomp;
	if($_ =~ m/INDEL/){
		next;
	}
	@m = split(/\s+/);
	@t = split(/;|=/,$m[7]);
	if($t[0] eq "DP"){
		$DP = $t[1];
	}
	else{
		$DP = $m[6];
	}
	$m[5] = int($DP/$m[6]*$m[5]);
	$m[6] = $DP;
	if($m[5]  <= $CUTOFF || $m[6] - $m[5] <= $CUTOFF){
		next;
	}

	if($m[6] > 1000){ ## repetitive region
		next;
	}
	#$m[6] = $DP;
	$#m = 6;
	print join("\t",@m),"\n";
}
