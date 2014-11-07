#!/usr/bin/perl
#N.tempfile  T.tempfile
$N = shift@ARGV;
open(I,"<$N");
while(<I>){
	chomp;
	@m = split(/\s+/);
	$hash{"$m[0]:$m[1]:$m[2]"}=$m[3];
	if($m[2] - $m[1] > 100){
		push @SUM, $m[3];
	}
}

@SUM_S = sort { $a <=> $b } @SUM;
$BASE = $SUM[int($#SUM_S*0.5)];
print $BASE,"\n";

$T = shift@ARGV;
open(I,"<$T");
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($hash{"$m[0]:$m[1]:$m[2]"} > 4){
		$m[3] = $m[3]/$hash{"$m[0]:$m[1]:$m[2]"}*$BASE;
	}
	print join("\t",@m),"\n";
}

