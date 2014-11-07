#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[2] - $m[1] < 200){
		next;
	}
	if($m[$#m] == 0){
		$S += $m[2] - $m[1];
		next;
	}
	if($m[3] + $m[4] != $m[8] + $m[9]){
		$IN += $m[$#m];
	}
}
print $S,"\t",$IN,"\n";

