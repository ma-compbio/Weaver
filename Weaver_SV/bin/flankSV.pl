#!/usr/bin/perl
#
$SV = shift@ARGV;
open(I,"<$SV");
while(<I>){
	chomp;
	@m = split(/\s+/);

	push @{$hash{$m[0]}}, $m[1];
	push @{$hash{$m[3]}}, $m[4];
=cut
	print $m[0],"\t",$m[1],"\t",$m[1]+10000,"\n";
	print $m[0],"\t",$m[1]-10000,"\t",$m[1],"\n";
	print $m[3],"\t",$m[4],"\t",$m[4]+10000,"\n";
	print $m[3],"\t",$m[4]-10000,"\t",$m[4],"\n";
=cut
}

foreach $chr (keys %hash){
	if($#{$hash{$chr}} > 0){
		for $i (1 .. $#{$hash{$chr}}){
			if(abs($hash{$chr}[$i] - $hash{$chr}[$i-1]) < 10000){
				$BLACK{$chr}{$hash{$chr}[$i]} = 1;
				$BLACK{$chr}{$hash{$chr}[$i-1]} = 1;
			}
		}
	}
}
open(S,">RETAINED");
open(I,"<$SV");
while(<I>){
	chomp;
	@m = split(/\s+/);
	if(!exists $BLACK{$m[0]}{$m[1]} && !exists $BLACK{$m[3]}{$m[4]}){
		print $m[0],"\t",$m[1],"\t",$m[1]+10000,"\n";
		print $m[0],"\t",$m[1]-10000,"\t",$m[1],"\n";
		print $m[3],"\t",$m[4],"\t",$m[4]+10000,"\n";
		print $m[3],"\t",$m[4]-10000,"\t",$m[4],"\n";
	}
}


