#!/usr/bin/perl
#
#1       109494847       -       3       110413400       +       120     0
open(I,"<",shift@ARGV);# reference sequence
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
print "##fileformat=VCFv4.2\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n";
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

