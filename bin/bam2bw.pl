#!/usr/bin/perl
use Parallel::ForkManager;
use FindBin qw($Bin);

$EXT = "$Bin/../external_bin";
#$pm=new Parallel::ForkManager();
$RUN_TYPE = shift@ARGV;
$bam = shift@ARGV;
$genome = "$bam.G1";
$genome_chr = "$bam.G2";
$P = shift@ARGV;
$SEX = shift@ARGV;

die "bam not found!\n" unless -e $bam;

system("$EXT/samtools view -H $bam | $Bin/bamHead2Genome.pl $SEX $bam");
open(I,"<$genome");
while(<I>){
	chomp;
	@m = split(/\s+/);
	push @ALL, $m[0];
	open(II,">$bam.$m[0].GG");
	print II $_,"\n";
}
$pm=new Parallel::ForkManager($P);
for($i = 0; $i < scalar @ALL; $i++){
	my $pid = $pm->start and next;
	$CHR = $ALL[$i];
	system("$EXT/samtools view -b $bam $CHR | $EXT/genomeCoverageBed -ibam stdin -g $bam.$CHR.GG -bga > $bam.$CHR.bed");
	system("cat $bam.$CHR.bed | $Bin/bedTofixWig.pl $bam.$CHR.GG > $bam.$CHR.wig; rm $bam.$CHR.bed");
	$pm->finish;
}
$pm->wait_all_children;
system("cat $bam.*.wig > $bam.wig; rm $bam.*.wig");
if($RUN_TYPE eq "lite"){
	system("rm $bam.*.GG");
}
else{
	system("$EXT/wigToBigWig $bam.wig $genome_chr $bam.bw; rm $bam.*.GG; rm $bam.G1 $bam.G2");
}

