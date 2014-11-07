#!/usr/bin/perl
#
#
$T_SNP = shift@ARGV;
open(I,"<",$T_SNP);
while(<I>){
	chomp;
	if(substr($_,0,1) eq "#"){
		next;
	}
	@m = split(/\s+/);
	$temp = "$m[0]:$m[1]:$m[3]:$m[4]";
	$hash{$temp}="$m[5]\t$m[6]";

}



$N_SNP = shift@ARGV;
open(I,"<",$N_SNP);
while(<I>){
	chomp;
	if(substr($_,0,1) eq "#"){
		next;
	}
	@m = split(/\s+/);
	$temp = "$m[0]:$m[1]:$m[3]:$m[4]";
	if(exists $hash{$temp}){
		print $_,"\t",$hash{$temp},"\n";
	}

}
