#!/usr/bin/perl
#
#./rmCEMT.pl REGION_CN_PHASE > xx
use FindBin qw($Bin);
$CENT = "$Bin/../data/CENT";
open(I,"<$CENT");
while(<I>){
	chomp;
	@m = split(/\s+/);
	$hash{$m[0]}[0] = $m[1];
	$hash{$m[0]}[1] = $m[2];
}
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if(substr($_,0,1) ne "c"){
		$CHR = "chr".$m[0];
	}
	else{
		$CHR = $m[0];
	}
	if($m[1] < $hash{$CHR}[0] && $m[2] <= $hash{$CHR}[1] && $m[2] > $hash{$CHR}[0]){
		print $m[0],"\t",$m[1],"\t",$hash{$CHR}[0],"\t",$m[3],"\t",$m[4],"\n";
	}
	elsif($m[2] >  $hash{$CHR}[1] && $m[1] < $hash{$CHR}[1] && $m[1] >= $hash{$CHR}[0]){
		print $m[0],"\t",$hash{$CHR}[1],"\t",$m[2],"\t",$m[3],"\t",$m[4],"\n";
	}
	elsif($m[1] < $hash{$CHR}[0] && $m[2] >  $hash{$CHR}[1]){
		print $m[0],"\t",$m[1],"\t",$hash{$CHR}[0],"\t",$m[3],"\t",$m[4],"\n";
		print $m[0],"\t",$hash{$CHR}[1],"\t",$m[2],"\t",$m[3],"\t",$m[4],"\n";
	}
	else{
		print $_,"\n";
	}
}

