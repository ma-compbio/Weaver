#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
#chr9    112319473       74      33      41      0       0       0
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[5] == 0 && $m[2] > 20 && $m[3] >= 6 && $m[4] >= 6){
		print $_,"\n";
	}

}

