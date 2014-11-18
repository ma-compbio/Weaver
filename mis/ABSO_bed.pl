#!/usr/bin/perl
#430625  1475054 TCGA-13-0727-01 8       3               331     19:24380818..28974637
#open(I,"<",shift@ARGV);
#grep ACOLD_p_TCGA_Batch17_SNP_N_GenomeWideSNP_6_D05_466158 ~/00686048/RUN/ABSOLUTE_CN | ~/00686048/deposit/Weaver/bin/ABSO_bed.pl > ABS
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	#@t = split(/:|\./,$m[$#m]);
	#print $t[0],"\t",$t[1],"\t",$t[3],"\t",$m[3],"\t",$m[4],"\n";
	print $m[1],"\t",$m[2],"\t",$m[3],"\t",$m[10],"\t",$m[11],"\n";
	print "\n";
}

