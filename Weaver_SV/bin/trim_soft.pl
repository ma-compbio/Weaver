#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
open(O1, ">SMALL_SOFT");
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[0] eq $m[3] && $m[2] ne $m[5] && abs($m[4] - $m[1]) < 4000){
		print O1 $_,"\n";
	}
	else{
		if($m[$#m] < 100){
			print $_,"\n";
		}
	}


}

