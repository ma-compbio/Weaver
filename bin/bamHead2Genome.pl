#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
#@SQ     SN:1    LN:249250621    UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta
$chr_flag = 0;
$number;
$SEX = shift@ARGV;
$BAM = shift@ARGV;

open(O1,">$BAM.G1");
open(O2,">$BAM.G2");
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[0] eq "\@SQ"){
		$chr = substr($m[1], 3);
		if($chr_flag == 0){
			$chr_flag = 1;
			if(substr($chr,0,1) eq "c"){
				$number = 0;
			}
			else{
				$number = 1;
			}
		}
		$length = substr($m[2], 3);
		if($number == 0){
			if(length($chr) > 5 || substr($chr,3,1) eq "M" || $SEX eq "F" && substr($chr,3,1) eq "Y"){
				next;
			}
			print O1 $chr,"\t",$length,"\n";
			print O2 $chr,"\t",$length,"\n";
		}
		else{
			if(length($chr) > 2 || substr($chr,0,1) eq "M" || $SEX eq "F" && substr($chr,0,1) eq "Y"){
				next;
			}
			print O1 $chr,"\t",$length,"\n";
			print O2 "chr".$chr,"\t",$length,"\n";
		}

	}
}

