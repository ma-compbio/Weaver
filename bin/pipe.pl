#!/usr/bin/perl
#
#$Bin = "/home/yangli9/tools/LUNA/HMM/LUNA";
use FindBin qw($Bin);
#$FLAG = shift@ARGV; # 0 for number 1 for chr
$EXT = "$Bin/../external_bin";
$RUN_TYPE = shift@ARGV;
$REGION = shift@ARGV;##GAP file with chr
$BAM = shift@ARGV;
$FA = shift@ARGV;#/home/yangli9/00686048/data/seq/hg19/number/noY/number_noY.fa
$P = shift@ARGV;# thread
#if($RUN_TYPE ne "lite"){
	$ONEKG = shift@ARGV;
#}
$SEX = shift@ARGV;

$TorN = shift@ARGV;

#system("samtools view -H $BAM | grep \"\\\@SQ\" | sed 's/^.*SN://g' | cut -f 1 |  xargs -I {} -n 1 -P 40 bash -c \"samtools mpileup -Q0 -d 100000 -uf $FA -r {} $BAM | bcftools view -vcg - > tmp.{}.vcf\"");
system("$EXT/samtools view -H $BAM | $Bin/bamHead2Genome_std.pl $SEX |  xargs -I {} -n 1 -P $P bash -c \"$EXT/samtools mpileup -Q0 -d 100000 -uf $FA -r {} $BAM | $EXT/bcftools view -vcg - > tmp.{}.vcf\"");
system("cat tmp*.vcf | $Bin/mergeVCF.pl > VCF; rm tmp*.vcf");
if($TorN eq "T" || !$TorN){
	system("$Bin/vcf2SNVlist.pl VCF | $Bin/clean_vcf.pl 8 | $Bin/check_GAP.pl $REGION > SNP");
}
if($TorN eq "N"){
	system("$Bin/vcf2SNVlist.pl VCF | $Bin/clean_vcf.pl 2 | $Bin/check_GAP.pl $REGION > SNP");
}
#$ONEKG/
system("zcat $Bin/../data/1000G_list | $Bin/corr.pl SNP > SNP_1000");# mv SNP_1000G SNP");
system("$Bin/density.pl SNP_1000 > SNP_dens");


use Parallel::ForkManager;
open(I,"<SNP_dens");
while(<I>){
	chomp;
	@m = split(/\s+/);
	$chr{$m[0]}=1;
}
@C = keys %chr;
if($P > scalar @C){
	$PP = scalar @C;
}
else{
	$PP = $P;
}


if($RUN_TYPE ne "lite"){
	$pm1=new Parallel::ForkManager($PP);
	for $t (0 .. $#C){
		my $pid1 = $pm1->start and next;
		system("$EXT/samtools view -h $BAM $C[$t] | $Bin/bamToSNPlink.pl | $Bin/newParse_memfix SNP_dens > LINK_$C[$t]");
		$pm1->finish;
	}
	$pm1->wait_all_children;
	system("cat LINK_* > LINK; rm LINK_*");

	$pm=new Parallel::ForkManager(2);
	for $i (0 .. 1){
		my $pid = $pm->start and next;
		if($i == 0){
			system("$Bin/SNPLinkCount.pl SNP_dens LINK > SNPLINK");
		}
		else{
			system("$Bin/ALL_1000G.pl SNP_dens $ONEKG $P");
		}
		$pm->finish;
	}
	$pm->wait_all_children;
}


