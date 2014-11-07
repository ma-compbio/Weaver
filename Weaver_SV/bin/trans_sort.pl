#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[0] gt $m[3]){
		print $m[3],"\t",$m[4],"\t",$m[5],"\t",$m[0],"\t",$m[1],"\t",$m[2],"\t",$m[6],"\n";
	}
	else{
		print $_,"\n";
	}

}
