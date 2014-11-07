#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[0] ne $m[3]  || ($m[0] eq $m[3] && $m[2] eq $m[5] && abs($m[1] - $m[4]) > 1000000) || ($m[0] eq $m[3] && $m[2] ne $m[5] && abs($m[1] - $m[4]) > 10000000)){
		print "chr",$m[0],"\t",$m[1],"\t",$m[2],"\tchr",$m[3],"\t",$m[4],"\t",$m[5],"\t",$m[6],"\t",$m[7],"\n";
		next;
	}
	#if($m[0] eq $m[3] && $m[2] eq $m[5])
}
