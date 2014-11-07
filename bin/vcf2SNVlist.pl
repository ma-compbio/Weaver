#!/usr/bin/perl
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	$LINE = $_;
	if(length($m[3]) > 1 || length($m[4]) > 1){
		next;
	}
	@s = split(/;/,$m[7]);
	foreach $key (@s){
		if(substr($key,0,3) eq "DP4"){
			$P = substr($key,4);
			@t = split(/,/,$P);
			if(!($t[0] == 0 && $t[1] == 0)){
				$m[5] = $t[2]+$t[3];## ALT ! not ref!
				$m[6] = $t[0]+$t[1]+$t[2]+$t[3];
				print join("\t",@m),"\n";
				last;
			}
		}
	}
}
