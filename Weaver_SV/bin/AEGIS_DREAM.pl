#!/usr/bin/perl
#

use FindBin qw($Bin);

#BEGIN {
#	$ENV{'LD_LIBRARY_PATH'}.="$Bin/../lib/";# loading bamtools lib
#	exec($^X, $0, @ARGV);
#}
$EXT = "$Bin/../../external_bin";
$BAM = shift@ARGV;
$BARCODE = $BAM.".AEGIS";
$FA_BWT = shift@ARGV; #just name
$FA = shift@ARGV;
$GAPFILE = "$Bin/../data/GAP_20140416";
#$GAPFILE = shift@ARGV; #"$Bin/../ManualGap";  ## with chr!!
$P = shift@ARGV;
use Parallel::ForkManager;
$pm=new Parallel::ForkManager(2);

=cut
if(!(-e $FA)){
	if(-e "$FA_BWT.fa"){
		$FA = "$FA_BWT.fa";
	}
	elsif(-e "$FA_BWT.fasta"){
		$FA = "$FA_BWT.fasta";
	}
	else{
		die "fasta reference not found! Provide by -F \n";
	}
}

die "$Bin/../../external_bin directory not working!" unless -x "$Bin/../../external_bin/bowtie";

die "GAP file not found!" unless -e $GAPFILE;

die "bam file not exist!" unless -e $BAM;


#system("export LD_LIBRARY_PATH=$Bin/../lib/:\$LD_LIBRARY_PATH");

for $i (0 .. 1){
	my $pid = $pm->start and next;
	if($i == 0){
		system("$Bin/Pair_bam $BAM > $BARCODE.NEAT");
		system("$Bin/superPair $BARCODE.NEAT > $BARCODE.PAIR");
	}
	else{
		system("$Bin/Bam_distri $BAM > $BARCODE.soft.fastq");
		system("$EXT/bowtie --suppress 6 --quiet --best -v 1 -p $P --un $BARCODE.soft.un.fastq $FA_BWT $BARCODE.soft.fastq $BARCODE.soft.bwt");
		system("$EXT/bwa mem -t $P $FA $BARCODE.soft.un.fastq > $BARCODE.soft.un.bwa.sam");
		system("cat $BARCODE.soft.un.bwa.sam | $Bin/format_bwa.pl > $BARCODE.bwa_partial_align");#change *
			system("cat $BARCODE.soft.bwt $BARCODE.bwa_partial_align | $Bin/breakpoint.pl > $BARCODE.all_soft");#change *
			system("sort -k 1,1 -k 2,2n $BARCODE.all_soft | $Bin/combine_soft.pl | sort -k 1,1 -k 2,2n | $Bin/combine_soft.pl | sort -k 1,1 -k 2,2n | $Bin/combine_soft.pl | sort -k 1,1 -k 2,2n > $BARCODE.SOFT");
		my $Num=0;
		my $N=0;
		while(1){
			open(I,"<$BARCODE.SOFT");
			while(<I>){
				$N++;
			}
			if($N == $Num){
				system("rm $BARCODE.SOFT2");
				last;
			}
			$Num = $N;
			$N=0;
			system("mv $BARCODE.SOFT $BARCODE.SOFT2");
			system("sort -k 1,1 -k 2,2n  $BARCODE.SOFT2 | $Bin/combine_soft.pl | sort -k 1,1 -k 2,2n > $BARCODE.SOFT");
			print $Num,"\t",$N,"\n";
		}
		system("$Bin/trans_sort.pl $BARCODE.SOFT > $BARCODE.SOFT_sort");

	}
	$pm->finish;
}
$pm->wait_all_children;
=cut

system("cat $BARCODE.SOFT_sort > $BARCODE.SOFT_sort_2");
system("$Bin/combineSuperPair $BARCODE.PAIR $BARCODE.SOFT_sort_2 > $BARCODE.FINAL_SV_");
system("cat $BARCODE.FINAL_SV_ > $BARCODE.FINAL_SV");
system("cat Only_SOFT > $BARCODE.FINAL_SOFT");
system("cat Only_Pair > $BARCODE.FINAL_PAIR"); ## 5 as cutoff
system("cat $BARCODE.FINAL_SV $BARCODE.FINAL_SOFT $BARCODE.FINAL_PAIR > $BARCODE.ALL_SV");## del or dup less than 2000 are discarted
