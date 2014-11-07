#!/usr/bin/perl
#
$K = shift@ARGV;
open(I,"<",$K);
while(<I>){
	chomp;
	@m = split(/\s+/);
	push @{$hash{$m[0]}},$m[1];
}

foreach $chr (keys %hash){
	for $i (0 .. 9){
		$s = 0;
		for $j ($i+1 .. $i+11){
			$s += abs(${$hash{$chr}}[$i] - ${$hash{$chr}}[$j]);
		}
		if($s/10 > 100000){
			$F{$chr}{${$hash{$chr}}[$i]}=1;
		}
	}
	for $i ($#{$hash{$chr}}-9 .. $#{$hash{$chr}}){
		$s = 0;
		for $j ($i-11 .. $i-1){
			$s += abs(${$hash{$chr}}[$i] - ${$hash{$chr}}[$j]);
		}
		if($s/10 > 100000){
			$F{$chr}{${$hash{$chr}}[$i]}=1;
		}
	}
	for $i (10 .. $#{$hash{$chr}}-10){
		$s1 = 0;
		$s2 = 0;
		for $j ($i-11 .. $i-1){
			$s1 += abs(${$hash{$chr}}[$i] - ${$hash{$chr}}[$j]);
		}
		for $j ($i+1 .. $i+11){
			$s2 += abs(${$hash{$chr}}[$i] - ${$hash{$chr}}[$j]);
		}
		#print $chr,"\t",${$hash{$chr}}[$i],"\t",$s/20,"\n";
		if($s1/10 > 100000 && $s2/10 > 100000){
			$F{$chr}{${$hash{$chr}}[$i]}=1;
		}
	}
}

open(I,"<",$K);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if(exists $F{$m[0]}{$m[1]}){
		#next;
		#print $_,"\n";
	}
	else{
		print $_,"\t0\n";
	}
}
