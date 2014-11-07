#!/usr/bin/perl
#

while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[4] < 40) {
		next;
	}
	if($m[5] =~ m/^(\d+)M$/){
		print $m[0],"\t",$m[2],"\t",$m[3],"\t",$1,"\t",$m[9],"\n";
	}
	if($m[5] =~ m/^(\d+)M(\d+)S$/){
		print $m[0],"\t",$m[2],"\t",$m[3],"\t",$1,"\t",$m[9],"\n";
	}
	if($m[5] =~ m/^(\d+)S(\d+)M$/){
		print $m[0],"\t",$m[2],"\t",$m[3],"\t",$2,"\t",substr($m[9],$1),"\n";
	}
}
