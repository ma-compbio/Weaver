#!/usr/bin/perl
#
#1       109494847       -       3       110413400       +       120     0
$REF = shift@ARGV;
open(I,"<",$REF);# reference sequence
while(<I>){
	chomp;
	if(substr($_,0,1) eq ">"){
		if($S ne ""){
			$hash{$chr} = $S;
			$S = "";
		}
		$chr = substr($_,1);
		next;
	}
	$S .= $_;
}
$hash{$chr} = $S;
$VERSION = 1;
print "##fileformat=VCFv4.2\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print "##source=AEGIS_V$VERSION\n##reference=file:$REF\n";
print "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n";
print "##FILTER=<ID=LowQual,Description=\"PE support below 3 or mapping quality below 20.\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">
##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">
##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">
##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">
##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

open(I,"<",shift@ARGV);## SV file
while(<I>){
	chomp;
	@m = split(/\s+/);
	$N++;
	$BASE1 = substr($hash{$m[0]}, $m[1], 1);
	$BASE2 = substr($hash{$m[3]}, $m[4], 1);
	if($m[2] eq "+" && $m[5] eq "+"){
		$M = $N+1;
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t$BASE1\]$m[3]:$m[4]\]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t$BASE2\]$m[0]:$m[1]\]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	if($m[2] eq "+" && $m[5] eq "-"){
		$M = $N+1;
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t$BASE1\[$m[3]:$m[4]\[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t\]$m[0]:$m[1]\]$BASE2\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	if($m[2] eq "-" && $m[5] eq "+"){
		$M = $N+1;
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t$BASE2\[$m[0]:$m[1]\[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t\]$m[3]:$m[4]\]$BASE1\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	if($m[2] eq "-" && $m[5] eq "-"){
		$M = $N+1;
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t\[$m[3]:$m[4]\[$BASE1\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t\[$m[0]:$m[1]\[$BASE2\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	$N = $M;
}

