#!/usr/bin/perl
#
#
#$KG = shift@ARGV; # 1 for chr 0 for no chr
$SNP = shift@ARGV;
open(I,"<",$SNP);
$S = <I>;
if(substr($S,0,1) eq "c"){
	$NUM_FLAG = 1;
}
else{
	$NUM_FLAG = 0;
}
#open(I,"<",$KG);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($NUM_FLAG == 1){
		$temp = "chr$m[0]:$m[1]:$m[2]:$m[3]";
	}
	if($NUM_FLAG == 0){
		$temp = "$m[0]:$m[1]:$m[2]:$m[3]";
	}
	$hash{$temp}=1;

}

open(I,"<",$SNP);
while(<I>){
	chomp;
	@m = split(/\s+/);
	$temp = "$m[0]:$m[1]:$m[3]:$m[4]";
	if(exists $hash{$temp}){
		print $_,"\n";
	}
}

