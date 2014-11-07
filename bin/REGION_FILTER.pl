#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[3]+$m[4] > 6 && abs($m[2] - $m[1]) < 10000){
		next;
	}
	else{
		print $_,"\n";
	}
}

