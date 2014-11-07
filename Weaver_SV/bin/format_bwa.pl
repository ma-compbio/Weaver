#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($#m < 5 ){
		next;
	}
	if($m[1] == 4){
		next;
	}
	if($m[1] == 0){
		$f = "+";
		print $m[0],"\t",$f,"\t",$m[2],"\t",$m[3],"\n";
	}
	if($m[1] == 16){
		$f = "-";
		print $m[0],"\t",$f,"\t",$m[2],"\t",$m[3],"\n";
	}



}
