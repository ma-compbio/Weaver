#!/usr/bin/perl
open(I,"<",shift@ARGV);## Tumor
while(<I>){
	chomp;
	if(substr($_,0,1) eq "#"){
		next;
	}
	@m = split(/\s+/);
	$LINE = $_;
	if(length($m[3]) > 1 || length($m[4]) > 1){
		next;
	}
	@s = split(/;/,$m[7]);
	foreach $key (@s){
		if(substr($key,0,3) eq "DP="){
			$DEPTH = substr($key,3);
		}
		if(substr($key,0,3) eq "DP4"){
			$P = substr($key,4);
			@t = split(/,/,$P);
			if(!($t[0] == 0 && $t[1] == 0)){
				$m[5] = $t[2]+$t[3];## ALT ! not ref!
					$m[6] = $t[0]+$t[1]+$t[2]+$t[3];
				last;
			}
		}
	}
	@{$hash{"$m[0]\t$m[1]"}} = ($DEPTH, $t[0]+$t[1],$t[2]+$t[3]);
}

open(I,"<",shift@ARGV);## Norm
while(<I>){
	chomp;
	if(substr($_,0,1) eq "#"){
		next;
	}
	@m = split(/\s+/);
	$LINE = $_;
	if(length($m[3]) > 1 || length($m[4]) > 1){
		next;
	}
	@s = split(/;/,$m[7]);
	foreach $key (@s){
		if(substr($key,0,3) eq "DP="){
			$DEPTH = substr($key,3);
		}
		if(substr($key,0,3) eq "DP4"){
			$P = substr($key,4);
			@t = split(/,/,$P);
			if(!($t[0] == 0 && $t[1] == 0)){
				$m[5] = $t[2]+$t[3];## ALT ! not ref!
					$m[6] = $t[0]+$t[1]+$t[2]+$t[3];
				last;
			}
		}
	}
	if(exists $hash{"$m[0]\t$m[1]"}){
		print "$m[0]\t$m[1]\t",join("\t",@{$hash{"$m[0]\t$m[1]"}}),"\t",$DEPTH,"\t",$t[0]+$t[1],"\t",$t[2]+$t[3],"\n";
		$S{"$m[0]\t$m[1]"}=1;
	}
	else{
		print "$m[0]\t$m[1]\t0\t0\t0\t",$DEPTH,"\t",$t[0]+$t[1],"\t",$t[2]+$t[3],"\n";
	}
}

foreach $key (keys %hash){
	if($S{$key} != 1){
		print "$key\t",join("\t",@{$hash{$key}}),"\t0\t0\t0\n";
	}
}
