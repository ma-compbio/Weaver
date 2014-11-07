#!/usr/bin/perl
#
use Parallel::ForkManager;

use FindBin qw($Bin);
=cut
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	push @ALL, $m[0];
	print $m[0],"\n";
}
=cut
for $i (1 .. 22){
	push @ALL, "chr$i";
}
push @ALL, "chrX";
#$FLAG = shift@ARGV;## 0 for number 1 for chr
$SNP = shift@ARGV;
$DIR_1000G = shift@ARGV;
$P = shift@ARGV;
open(I,"<",$SNP);
$S = <I>;
if(substr($S,0,1) eq "c"){
	$FLAG = 1;
}
else{
	$FLAG = 0;
}

$pm=new Parallel::ForkManager($P);
for $i (0 .. $#ALL){
	my $pid = $pm->start and next;
	print $i,"\n";
	$chr = $ALL[$i];
	print $chr,"\n";
	#system("proce.pl ALL.$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf > $chr.1000G");
	system("$Bin/LD.pl $FLAG  $chr $SNP $DIR_1000G/$chr.1000G > $chr.1000G.link");
	$pm->finish;
}
$pm->wait_all_children;
system("cat chr*.1000G.link > SNP_LINK_1000G; rm *.1000G.link");


