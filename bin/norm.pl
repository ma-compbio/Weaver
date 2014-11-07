#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	push @ALL, $m[$#m];
	$hash{"$m[0]:$m[1]:$m[2]"} = $m[3];
}
@N = sort{$a<=>$b} @ALL;
$F = $N[int($#N/2)];


open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	#push @ALL, $m[$#m];
	if($hash{"$m[0]:$m[1]:$m[2]"} > 5){
		$m[3] = $m[3]/$hash{"$m[0]:$m[1]:$m[2]"}*$F;
	}
	print join("\t",@m),"\n";
}



