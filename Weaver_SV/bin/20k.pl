#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
$DIS = shift@ARGV;
$PAIR_LIMIT = shift@ARGV;
$SOFT_LIMIT = shift@ARGV;
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[0] ne $m[3] || $m[2] eq $m[5] || ($m[0] eq $m[3] && abs($m[1]-$m[4]) > $DIS)){
		if($m[7] == 0 && $m[6] < $PAIR_LIMIT){
			next;
		}
		if($m[6] == 0 && $m[7] < $SOFT_LIMIT){
			next;
		}
		print $_,"\n";
	}
}

